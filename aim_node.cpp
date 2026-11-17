// RoboMaster Aim Lab 学员自瞄节点
//
// 处理链路：订阅装甲板角点/云台状态/射击参数 → 选正对相机的板 → 轨迹拟合与运动预测
//          → 抛体解算求 yaw/pitch → 到位开火（含收尾控火）→ 发布控制指令。
// 支持 4/6/10 m 固定靶与 6 m 平移、旋转、平移+旋转组合靶。
#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include <rclcpp/rclcpp.hpp>
#include <rm_aim_sim_bridge/msg/aim_command.hpp>
#include <rm_aim_sim_bridge/msg/armor_corners_array.hpp>
#include <rm_aim_sim_bridge/msg/gimbal_state.hpp>
#include <rm_aim_sim_bridge/msg/telemetry.hpp>

using namespace std::chrono_literals;
using ArmorsMsg = rm_aim_sim_bridge::msg::ArmorCornersArray;
using GimbalMsg = rm_aim_sim_bridge::msg::GimbalState;
using TelemetryMsg = rm_aim_sim_bridge::msg::Telemetry;
using CommandMsg = rm_aim_sim_bridge::msg::AimCommand;

namespace {

struct Vec3 {
  double x = 0.0, y = 0.0, z = 0.0;
  Vec3 operator+(const Vec3 & o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3 & o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
  double norm() const { return std::hypot(x, y, z); }
};

double wrapPi(double a) {
  while (a > M_PI) a -= 2.0 * M_PI;
  while (a < -M_PI) a += 2.0 * M_PI;
  return a;
}

double stampToSec(const builtin_interfaces::msg::Time & t) {
  return static_cast<double>(t.sec) + static_cast<double>(t.nanosec) * 1e-9;
}

struct Sample { double t; Vec3 p; };

bool solve3(const double A[9], const double b[3], double x[3]);

// 3x3 高斯消元（列主元），供二次拟合与圆拟合共用。
bool solve3(const double A[9], const double b[3], double x[3]) {
  double a[3][3] = {{A[0], A[1], A[2]}, {A[3], A[4], A[5]}, {A[6], A[7], A[8]}};
  double r[3] = {b[0], b[1], b[2]};
  for (int col = 0; col < 3; ++col) {
    int piv = col;
    for (int i = col + 1; i < 3; ++i)
      if (std::fabs(a[i][col]) > std::fabs(a[piv][col])) piv = i;
    if (std::fabs(a[piv][col]) < 1e-12) return false;
    if (piv != col) {
      for (int c = 0; c < 3; ++c) std::swap(a[piv][c], a[col][c]);
      std::swap(r[piv], r[col]);
    }
    for (int i = col + 1; i < 3; ++i) {
      const double f = a[i][col] / a[col][col];
      for (int c = col; c < 3; ++c) a[i][c] -= f * a[col][c];
      r[i] -= f * r[col];
    }
  }
  for (int i = 2; i >= 0; --i) {
    double s = r[i];
    for (int c = i + 1; c < 3; ++c) s -= a[i][c] * x[c];
    x[i] = s / a[i][i];
  }
  return true;
}

// 梯形速度 1D 数值预测：平移靶做「加速 → 巡航 → 减速 → 换向」往返运动
// （speed=1.5 m/s、acceleration=3 m/s²、travel=1.5 m），分段精确外推，换向点不滞后。
double predictTrapezoid1D(double y, double v, double lead,
                          double peak, double accel, double travel) {
  double yy = y, vv = v, t = 0.0;
  const double dt = 0.004;
  while (t < lead - 1e-12) {
    const double step = std::min(dt, lead - t);
    double a = 0.0;
    const double stop = vv * vv / (2.0 * accel);  // 减速到 0 所需距离
    if (std::fabs(vv) < 1e-3) {
      if (yy >= travel - 1e-3) a = -accel;         // + 边界换向
      else if (yy <= -travel + 1e-3) a = accel;    // - 边界换向
    } else if (vv > 0.0) {
      if (travel - yy <= stop) a = -accel;         // 减速换向
      else if (vv < peak) a = accel;               // 加速
    } else {
      if (yy + travel <= stop) a = accel;          // 减速换向
      else if (vv > -peak) a = -accel;             // 加速
    }
    yy += vv * step + 0.5 * a * step * step;
    vv += a * step;
    yy = std::clamp(yy, -travel, travel);
    t += step;
  }
  return yy;
}

// Kåsa 最小二乘圆拟合：由轨迹样本估计旋转中心与半径（供旋转/组合靶外推）。
bool fitCircle(const std::deque<Sample> & smp, double & cx, double & cy, double & r) {
  if (smp.size() < 4) return false;
  double Sx = 0, Sy = 0, Sx2 = 0, Sy2 = 0, Sxy = 0, Sz = 0, Sxz = 0, Syz = 0;
  for (const auto & s : smp) {
    const double x = s.p.x, y = s.p.y, z = x * x + y * y;
    Sx += x; Sy += y; Sx2 += x * x; Sy2 += y * y; Sxy += x * y;
    Sz += z; Sxz += x * z; Syz += y * z;
  }
  const double A[9] = {Sx2, Sxy, Sx, Sxy, Sy2, Sy, Sx, Sy, static_cast<double>(smp.size())};
  const double b[3] = {-Sxz, -Syz, -Sz};
  double x[3] = {0, 0, 0};
  if (!solve3(A, b, x)) return false;
  cx = -x[0] / 2.0;
  cy = -x[1] / 2.0;
  r = std::sqrt(std::max(0.0, cx * cx + cy * cy - x[2]));
  return std::isfinite(cx) && std::isfinite(cy) && std::isfinite(r);
}

// 最小二乘二次拟合：在最新样本时刻求 p = a + b·τ + c·τ² 的位置/速度/加速度。
bool fitQuadratic(const std::deque<Sample> & smp, Vec3 & pos, Vec3 & vel, Vec3 & acc) {
  if (smp.size() < 4) return false;
  const double tmid = (smp.front().t + smp.back().t) * 0.5;
  double S[5] = {0, 0, 0, 0, 0};
  double P[3] = {0, 0, 0}, Tp[3] = {0, 0, 0}, T2p[3] = {0, 0, 0};
  for (const auto & s : smp) {
    const double tau = s.t - tmid;
    const double t2 = tau * tau;
    S[0] += 1; S[1] += tau; S[2] += t2; S[3] += t2 * tau; S[4] += t2 * t2;
    const double px[3] = {s.p.x, s.p.y, s.p.z};
    for (int k = 0; k < 3; ++k) {
      P[k] += px[k]; Tp[k] += tau * px[k]; T2p[k] += t2 * px[k];
    }
  }
  const double A[9] = {S[0], S[1], S[2], S[1], S[2], S[3], S[2], S[3], S[4]};
  double a[3] = {0, 0, 0}, b[3] = {0, 0, 0}, c[3] = {0, 0, 0};
  for (int k = 0; k < 3; ++k) {
    const double rhs[3] = {P[k], Tp[k], T2p[k]};
    double x[3] = {0, 0, 0};
    if (!solve3(A, rhs, x)) return false;
    a[k] = x[0]; b[k] = x[1]; c[k] = x[2];
  }
  // 在最新样本时刻求值（中点求值会引入约半个窗口的滞后）。
  const double tl = smp.back().t - tmid;
  for (int k = 0; k < 3; ++k) {
    const double p = a[k] + b[k] * tl + c[k] * tl * tl;
    const double v = b[k] + 2.0 * c[k] * tl;
    const double acc_k = 2.0 * c[k];
    if (k == 0) { pos.x = p; vel.x = v; acc.x = acc_k; }
    else if (k == 1) { pos.y = p; vel.y = v; acc.y = acc_k; }
    else { pos.z = p; vel.z = v; acc.z = acc_k; }
  }
  return true;
}

}  // namespace

class AimNode : public rclcpp::Node {
 public:
  AimNode() : Node("aim_node") {
    yaw_tol_scale_ = declare_parameter<double>("yaw_tol_scale", 0.7);
    pitch_tol_scale_ = declare_parameter<double>("pitch_tol_scale", 0.7);
    lead_offset_ = declare_parameter<double>("lead_offset", 0.03);

    const auto qos = rclcpp::QoS(rclcpp::KeepLast(1)).reliable();
    armors_sub_ = create_subscription<ArmorsMsg>(
        "/aim_sim/armors", qos,
        [this](ArmorsMsg::ConstSharedPtr m) {
          std::lock_guard<std::mutex> l(mu_);
          armors_ = m;
        });
    gimbal_sub_ = create_subscription<GimbalMsg>(
        "/aim_sim/gimbal", qos,
        [this](GimbalMsg::ConstSharedPtr m) {
          std::lock_guard<std::mutex> l(mu_);
          gimbal_ = m;
        });
    telemetry_sub_ = create_subscription<TelemetryMsg>(
        "/aim_sim/telemetry", qos,
        [this](TelemetryMsg::ConstSharedPtr m) {
          std::lock_guard<std::mutex> l(mu_);
          telemetry_ = m;
        });
    command_pub_ = create_publisher<CommandMsg>("/aim_sim/command", qos);
    timer_ = create_wall_timer(20ms, [this] { onTimer(); });

    RCLCPP_INFO(get_logger(), "aim_node started (quadratic fit + circle predict + endgame control)");
  }

 private:
  // ---- 主循环 ----
  void onTimer() {
    std::lock_guard<std::mutex> lock(mu_);
    if (!armors_) return;  // 尚未收到观测，无法构造同源指令

    CommandMsg cmd;
    cmd.header = armors_->header;
    cmd.session_id = armors_->session_id;

    if (armors_->session_id != track_session_) {
      track_session_ = armors_->session_id;
      resetTracking();
    }

    const auto & armors = armors_->armors;
    const double v = telemetry_ ? telemetry_->speed : 23.0;
    const double g = telemetry_ ? telemetry_->gravity : 9.7833;
    const double L = telemetry_ ? telemetry_->barrel_length : 0.22;
    const double delay = telemetry_ ? telemetry_->delay : 0.015;
    const double now = stampToSec(armors_->header.stamp);

    if (armors.empty()) {
      handleNoTarget(now, cmd);
      command_pub_->publish(cmd);
      return;
    }
    settleLightOff(now);  // 目标重现：结算本次空窗是否为真命中灭灯

    const Vec3 center = selectArmor(armors);
    pushSample(center, now);
    const Vec3 aim_point = predict(center, v, delay);

    double yaw = 0.0, pitch = 0.0;
    solveAim(aim_point, v, g, L, yaw, pitch);

    last_seen_t_ = now;
    last_yaw_ = yaw;
    last_pitch_ = pitch;

    cmd.yaw = yaw;
    cmd.pitch = pitch;
    cmd.enable = true;
    cmd.fire = fireGate(aim_point, yaw, pitch);
    applyEndgame(now, cmd);           // 收尾控火：减少击毁后在途弹丸的计分浪费
    recordShot(now, cmd.fire, center, v);  // 记录发射时刻，供命中灭灯判定

    command_pub_->publish(cmd);

    if ((get_clock()->now() - last_log_).seconds() > 1.0) {
      last_log_ = get_clock()->now();
      RCLCPP_INFO(get_logger(),
          "c=(%.2f,%.2f) w=%.2f hits=%d tgt=(%.4f,%.4f) fire=%d fit=%d",
          center.x, center.y, omega_, hit_count_, yaw, pitch,
          cmd.fire, fit_valid_ ? 1 : 0);
    }
  }

  void resetTracking() {
    samples_.clear();
    fit_valid_ = false;
    pos_ = vel_ = acc_ = {0, 0, 0};
    omega_ = 0.0;
    vang_valid_ = false;
    center_valid_ = false;
    circle_r_ = 0.0;
    center_vel_ = {0, 0, 0};
    cc_valid_ = false;
    last_seen_t_ = -1e9;
    hit_count_ = 0;
    light_off_ = false;
    last_burst_t_ = -1e9;
    endgame_off_ = false;
    shot_times_.clear();
    was_firing_ = false;
  }

  // ---- 无目标：短空窗（命中灭灯/渲染丢帧）保持上一瞄准角，长空窗才搜索 ----
  void handleNoTarget(double now, CommandMsg & cmd) {
    if (now - last_seen_t_ >= kSearchDelay) {
      searchMode(cmd);
      return;
    }
    if (!light_off_) {
      light_off_ = true;
      light_off_start_ = now;
    }
    cmd.yaw = last_yaw_;
    cmd.pitch = last_pitch_;
    cmd.enable = true;
    cmd.fire = false;
  }

  // 灭灯结算：真命中灭灯发生在「发射 + 飞行时间」附近；渲染丢帧空窗随机，不在窗口内。
  // 该判定在 armors 恢复时（onTimer 的目标分支开头）调用。
  void settleLightOff(double now) {
    if (!light_off_) return;
    light_off_ = false;
    const double dur = now - light_off_start_;
    if (dur < 0.02 || dur > 0.15) return;
    for (const double st : shot_times_) {
      if (std::fabs((now - st) - last_t_flight_) < 0.12) {
        ++hit_count_;
        break;
      }
    }
  }

  // ---- 选板：四角叉积求朝外法向，取法向 x 最小（-x 面，正对云台）----
  Vec3 selectArmor(const std::vector<rm_aim_sim_bridge::msg::ArmorCorners> & armors) const {
    Vec3 best{0, 0, 0};
    double best_nx = 1e9;
    for (const auto & a : armors) {
      Vec3 c, n;
      plateGeometry(a, c, n);
      if (n.x < best_nx) { best_nx = n.x; best = c; }
    }
    return best;
  }

  static void plateGeometry(const rm_aim_sim_bridge::msg::ArmorCorners & a, Vec3 & c, Vec3 & n) {
    c = {0, 0, 0};
    for (const auto & p : a.corners) { c.x += p.x; c.y += p.y; c.z += p.z; }
    c = c / 4.0;
    const auto & tl = a.corners[0];
    const auto & tr = a.corners[1];
    const auto & bl = a.corners[3];
    const Vec3 v1{tr.x - tl.x, tr.y - tl.y, tr.z - tl.z};
    const Vec3 v2{bl.x - tl.x, bl.y - tl.y, bl.z - tl.z};
    n = {v2.y * v1.z - v2.z * v1.y,
         v2.z * v1.x - v2.x * v1.z,
         v2.x * v1.y - v2.y * v1.x};
    const double len = n.norm();
    if (len > 1e-9) n = n / len;
  }

  // ---- 轨迹更新：窗口维护 + 圆拟合 + 二次拟合 + 角速度估计 ----
  void pushSample(Vec3 c, double now) {
    // 观测跳变：旋转靶换板清空窗口；平移靶正对板被挡则用上一帧状态桥接保持连续。
    if (!samples_.empty()) {
      const double d = (c - samples_.back().p).norm();
      if (d > kJumpThreshold || now <= samples_.back().t) {
        if (std::fabs(omega_) > 0.6) {
          samples_.clear();
          fit_valid_ = false;
        } else if (fit_valid_) {
          c = pos_ + vel_ * (now - samples_.back().t);
        }
      }
    }
    samples_.push_back({now, c});
    while (!samples_.empty() && now - samples_.front().t > kWindowSec) samples_.pop_front();
    if (samples_.size() > kMaxSamples) samples_.pop_front();

    // 圆拟合：旋转中心与半径跨换板保持，圆心平移速度供组合靶外推。
    double ccx = 0, ccy = 0, cr = 0;
    if (fitCircle(samples_, ccx, ccy, cr) && cr < 3.0) {
      if (cc_valid_ && now > last_cc_t_) {
        const double dtc = now - last_cc_t_;
        if (dtc > 1e-3 && dtc < 0.5) {
          const double vcx = (ccx - last_ccx_) / dtc;
          const double vcy = (ccy - last_ccy_) / dtc;
          if (std::hypot(vcx, vcy) < 5.0) {
            center_vel_.x = center_vel_.x * 0.6 + vcx * 0.4;
            center_vel_.y = center_vel_.y * 0.6 + vcy * 0.4;
          }
        }
      }
      last_ccx_ = ccx; last_ccy_ = ccy; last_cc_t_ = now; cc_valid_ = true;
      if (center_valid_) {
        center_.x = center_.x * 0.7 + ccx * 0.3;
        center_.y = center_.y * 0.7 + ccy * 0.3;
        circle_r_ = circle_r_ * 0.7 + cr * 0.3;
      } else {
        center_.x = ccx; center_.y = ccy; circle_r_ = cr;
        center_valid_ = true;
      }
    }

    fit_valid_ = fitQuadratic(samples_, pos_, vel_, acc_);
    if (fit_valid_) {
      // 角速度 = 速度方向变化率（平移靶方向恒定 → ω≈0，比 (v×a)/|v|² 抗噪）。
      const double vsp = std::hypot(vel_.x, vel_.y);
      if (vsp > 0.3) {
        const double th = std::atan2(vel_.y, vel_.x);
        if (vang_valid_ && now > last_vang_t_) {
          const double dt = now - last_vang_t_;
          const double w_inst = wrapPi(th - vang_) / dt;
          // 旋转方向恒定：已建立 |ω|>1.5 后忽略反向跳变（换板时速度方向突变 90°）。
          const bool sign_flip = std::fabs(omega_) > 1.5 && (omega_ * w_inst) < 0.0;
          if (dt > 1e-3 && dt < 0.5 && std::fabs(w_inst) < 15.0 && !sign_flip)
            omega_ = omega_ * 0.6 + w_inst * 0.4;
        }
        vang_ = th;
        last_vang_t_ = now;
        vang_valid_ = true;
      } else {
        vang_valid_ = false;
      }
    } else {
      vel_ = acc_ = {0, 0, 0};  // omega_ 与圆参数保持
    }
  }

  // ---- 运动预测 ----
  Vec3 predict(const Vec3 & c, double v, double delay) const {
    const double lead = c.norm() / v + delay + lead_offset_;
    if (!fit_valid_) return c;
    // 旋转/组合靶：绕圆心外推，圆心随平移速度推进。
    if (std::fabs(omega_) > 0.6 && center_valid_) {
      const double ccx = center_.x + center_vel_.x * lead;
      const double ccy = center_.y + center_vel_.y * lead;
      const double th = omega_ * lead;
      const double rx = c.x - ccx, ry = c.y - ccy;
      const double cs = std::cos(th), sn = std::sin(th);
      return {ccx + rx * cs - ry * sn,
              ccy + rx * sn + ry * cs,
              pos_.z + vel_.z * lead};
    }
    // 平移/固定靶：y 梯形速度模型，x/z 二阶外推。
    return {pos_.x + vel_.x * lead + 0.5 * acc_.x * lead * lead,
            predictTrapezoid1D(pos_.y, vel_.y, lead, kPeakSpeed, kAccel, kTravel),
            pos_.z + vel_.z * lead};
  }

  // ---- 抛体解算：从枪口（云台前方 L）打到目标，迭代求 yaw/pitch（低伸弹道）----
  void solveAim(const Vec3 & c, double v, double g, double L,
                double & yaw, double & pitch) const {
    yaw = std::atan2(c.y, c.x);
    pitch = std::atan2(c.z, std::hypot(c.x, c.y));
    for (int it = 0; it < 4; ++it) {
      const double cp = std::cos(pitch), sp = std::sin(pitch);
      const Vec3 muzzle{L * cp * std::cos(yaw), L * cp * std::sin(yaw), L * sp};
      const Vec3 d = c - muzzle;
      yaw = std::atan2(d.y, d.x);
      const double h = std::hypot(d.x, d.y);
      if (h < 1e-6) { pitch = std::atan2(d.z, h); continue; }
      const double v2 = v * v;
      const double disc = v2 * v2 - g * (g * h * h + 2.0 * d.z * v2);
      pitch = (disc < 0.0) ? std::atan2(d.z, h)
                           : std::atan((v2 - std::sqrt(disc)) / (g * h));
    }
  }

  // ---- 开火门限 ----
  bool fireGate(const Vec3 & c, double target_yaw, double target_pitch) const {
    if (!fit_valid_) return false;
    // 旋转靶：预测精度受换板影响，保留云台反馈门限过滤不准时刻。
    if (std::fabs(omega_) > 0.6 && gimbal_) {
      const double dist = c.norm();
      const double yaw_tol = std::atan2(kHalfWidth, dist) * yaw_tol_scale_;
      const double pitch_tol = std::atan2(kHalfHeight, dist) * pitch_tol_scale_;
      const double ang_rate = std::hypot(vel_.x, vel_.y) / dist;
      const double lag_tol = ang_rate * kGateLag;
      return std::fabs(wrapPi(target_yaw - gimbal_->yaw)) <= yaw_tol + lag_tol &&
             std::fabs(target_pitch - gimbal_->pitch) <= pitch_tol + lag_tol;
    }
    // 平移/固定靶：仿真为理想云台（指令即实际角），反馈比较只会误关换向点的火。
    return true;
  }

  // ---- 收尾控火：命中数达阈值后停火让在途弹丸完成击毁，不足则脉冲补发 ----
  void applyEndgame(double now, CommandMsg & cmd) {
    if (hit_count_ < kStopHits || endgame_off_) return;
    if (last_burst_t_ < -1e8) last_burst_t_ = now;
    const double since = now - last_burst_t_;
    if (since < 0.5)
      cmd.fire = false;  // 静默期：等在途弹丸命中
    else
      cmd.fire = cmd.fire && (std::fmod(since - 0.5, 0.5) < 0.2);  // 补发脉冲
    if (since > 2.0) endgame_off_ = true;  // 低命中率场景恢复连续开火
  }

  // 估计发射时刻队列（fire 上升沿 1 发 + 持续 fire 每 0.09 s 一发）。
  void recordShot(double now, bool fire, const Vec3 & c, double v) {
    if (fire && !was_firing_) {
      shot_times_.push_back(now);
      last_shot_est_t_ = now;
    } else if (fire && now - last_shot_est_t_ >= 0.09) {
      shot_times_.push_back(now);
      last_shot_est_t_ = now;
    }
    was_firing_ = fire;
    last_t_flight_ = c.norm() / v;
    while (!shot_times_.empty() && now - shot_times_.front() > last_t_flight_ + 0.4)
      shot_times_.pop_front();
  }

  void searchMode(CommandMsg & cmd) const {
    const double t = std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
    cmd.yaw = 0.9 * std::sin(1.2 * t);
    cmd.pitch = 0.0;
    cmd.enable = true;
    cmd.fire = false;
  }

  // ---- 常量 ----
  static constexpr double kHalfWidth = 0.0675;    // 装甲板半宽（m）
  static constexpr double kHalfHeight = 0.0625;   // 装甲板半高（m）
  static constexpr int kStopHits = 13;            // 收尾停火命中数 ≈ 18 − 在途弹丸 − 检测延迟
  static constexpr double kGateLag = 0.06;        // 旋转靶门限的角速度×延迟补偿（s）
  static constexpr double kWindowSec = 0.3;       // 轨迹样本时间窗（s）
  static constexpr int kMaxSamples = 12;          // 轨迹样本数上限
  static constexpr double kJumpThreshold = 0.35;  // 观测跳变阈值（m）
  static constexpr double kSearchDelay = 0.3;     // 丢失多久转入搜索（s）
  static constexpr double kPeakSpeed = 1.5;       // 平移靶峰值速度（m/s）
  static constexpr double kAccel = 3.0;           // 平移靶加速度（m/s²）
  static constexpr double kTravel = 1.5;          // 平移靶行程（m）

  // ---- 参数 ----
  double yaw_tol_scale_{0.7};
  double pitch_tol_scale_{0.7};
  double lead_offset_{0.03};

  // ---- 轨迹状态 ----
  std::string track_session_;
  std::deque<Sample> samples_;
  bool fit_valid_{false};
  Vec3 pos_{0, 0, 0}, vel_{0, 0, 0}, acc_{0, 0, 0};
  double omega_{0.0};
  double vang_{0.0};
  bool vang_valid_{false};
  double last_vang_t_{0.0};
  Vec3 center_{0, 0, 0};           // 旋转中心
  bool center_valid_{false};
  double circle_r_{0.0};
  Vec3 center_vel_{0, 0, 0};       // 旋转中心平移速度（组合靶）
  bool cc_valid_{false};
  double last_ccx_{0.0}, last_ccy_{0.0}, last_cc_t_{0.0};

  // ---- 目标丢失 / 灭灯计数状态 ----
  double last_seen_t_{-1e9};
  double last_yaw_{0.0}, last_pitch_{0.0};
  bool light_off_{false};
  double light_off_start_{0.0};
  int hit_count_{0};
  bool endgame_off_{false};
  double last_burst_t_{-1e9};
  std::deque<double> shot_times_;
  double last_shot_est_t_{-1e9};
  bool was_firing_{false};
  double last_t_flight_{0.2};

  // ---- ROS 通信 ----
  std::mutex mu_;
  ArmorsMsg::ConstSharedPtr armors_{nullptr};
  GimbalMsg::ConstSharedPtr gimbal_{nullptr};
  TelemetryMsg::ConstSharedPtr telemetry_{nullptr};
  rclcpp::Subscription<ArmorsMsg>::SharedPtr armors_sub_;
  rclcpp::Subscription<GimbalMsg>::SharedPtr gimbal_sub_;
  rclcpp::Subscription<TelemetryMsg>::SharedPtr telemetry_sub_;
  rclcpp::Publisher<CommandMsg>::SharedPtr command_pub_;
  rclcpp::TimerBase::SharedPtr timer_;
  rclcpp::Time last_log_{0, 0, RCL_ROS_TIME};
};

int main(int argc, char ** argv) {
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<AimNode>());
  rclcpp::shutdown();
  return 0;
}
