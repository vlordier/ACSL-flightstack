#include "flightstack/ardupilot_interface/ardupilot_actuator.hpp"
#include <mavros_msgs/msg/actuator_control.hpp>
#include <mavros_msgs/srv/command_bool.hpp>
#include <mavros_msgs/srv/set_mode.hpp>
#include <mavros_msgs/msg/state.hpp>

void ArduPilotActuator::arm()
{
    auto client = node_->create_client<mavros_msgs::srv::CommandBool>(
        ns_ + "/mavros/cmd/arming");
    auto request = std::make_shared<mavros_msgs::srv::CommandBool::Request>();
    request->value = true;
    client->async_send_request(request);
    RCLCPP_INFO(node_->get_logger(), "[ardupilot] Arm");
}

void ArduPilotActuator::disarm()
{
    auto client = node_->create_client<mavros_msgs::srv::CommandBool>(
        ns_ + "/mavros/cmd/arming");
    auto request = std::make_shared<mavros_msgs::srv::CommandBool::Request>();
    request->value = false;
    client->async_send_request(request);
    RCLCPP_INFO(node_->get_logger(), "[ardupilot] Disarm");
}

void ArduPilotActuator::set_guided_mode()
{
    auto client = node_->create_client<mavros_msgs::srv::SetMode>(
        ns_ + "/mavros/set_mode");
    auto request = std::make_shared<mavros_msgs::srv::SetMode::Request>();
    request->custom_mode = "GUIDED";
    client->async_send_request(request);
    RCLCPP_INFO(node_->get_logger(), "[ardupilot] GUIDED mode");
}

void ArduPilotActuator::publish_actuator_controls(const std::array<double, 4>& controls)
{
    auto msg = mavros_msgs::msg::ActuatorControl();
    msg.group_mix = 0;
    for (size_t i = 0; i < 4 && i < controls.size(); ++i)
        msg.controls[i] = controls[i];
    actuator_pub_->publish(msg);
}
