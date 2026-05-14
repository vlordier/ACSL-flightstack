#ifndef FLIGHTSTACK__ARDUPLIOT_ACTUATOR_HPP_
#define FLIGHTSTACK__ARDUPLIOT_ACTUATOR_HPP_

#include <rclcpp/rclcpp.hpp>
#include <mavros_msgs/msg/actuator_control.hpp>
#include <array>
#include <string>

class ArduPilotActuator {
public:
    ArduPilotActuator(rclcpp::Node::SharedPtr node, const std::string& ns)
        : node_(node), ns_(ns)
    {
        actuator_pub_ = node_->create_publisher<mavros_msgs::msg::ActuatorControl>(
            ns_ + "/mavros/actuator_control", 10);
        state_sub_ = node_->create_subscription<mavros_msgs::msg::State>(
            ns_ + "/mavros/state", 10,
            [this](const mavros_msgs::msg::State::SharedPtr msg) {
                current_state_ = *msg;
            });
    }

    void arm();
    void disarm();
    void set_guided_mode();
    void publish_actuator_controls(const std::array<double, 4>& controls);
    bool is_armed() const { return current_state_.armed; }

private:
    rclcpp::Node::SharedPtr node_;
    std::string ns_;
    rclcpp::Publisher<mavros_msgs::msg::ActuatorControl>::SharedPtr actuator_pub_;
    rclcpp::Subscription<mavros_msgs::msg::State>::SharedPtr state_sub_;
    mavros_msgs::msg::State current_state_;
};

#endif
