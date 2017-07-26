#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *

print "Initializing..."
# Function to create a DH Transform Matrix and substitute symbol values from dictionary
def dhTransformMatrix(alpha, a, d, q, s):
    TF = Matrix([[            cos(q),           -sin(q),           0,             a],
                 [ sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                 [ sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                 [                  0,                0,           0,             1]])
    return TF.subs(s)

print "Creating symbols and dictionary..."
# Define DH param symbols
q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') #theta_i
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

# Dictionary with DH param values
s = {alpha0:     0,  a0:      0, d1: 0.75, 
     alpha1: -pi/2,  a1:   0.35, d2:    0,  q2: q2-pi/2,  
     alpha2:     0,  a2:   1.25, d3:    0,
     alpha3: -pi/2,  a3:  0.054, d4: 1.50,
     alpha4:  pi/2,  a4:      0, d5:    0,
     alpha5: -pi/2,  a5:      0, d6:    0,
     alpha6:     0,  a6:      0, d7: 0.2305, q7: 0} #0.3405 0.303 0.2305 0.193

print "Creating individual transformation matrices..."
# Create individual transformation matrices
T0_1 = dhTransformMatrix(alpha0, a0, d1, q1, s)
T1_2 = dhTransformMatrix(alpha1, a1, d2, q2, s)
T2_3 = dhTransformMatrix(alpha2, a2, d3, q3, s)
T3_4 = dhTransformMatrix(alpha3, a3, d4, q4, s)
T4_5 = dhTransformMatrix(alpha4, a4, d5, q5, s)
T5_6 = dhTransformMatrix(alpha5, a5, d6, q6, s)
T6_G = dhTransformMatrix(alpha6, a6, d7, q7, s)

print "Building transformation matrix T0_G..."
# Build total transformation matrix by post-multiplying individual transforms
T0_2 = simplify(T0_1 * T1_2)
T0_3 = simplify(T0_2 * T2_3)
T0_4 = simplify(T0_3 * T3_4)
T0_5 = simplify(T0_4 * T4_5)
T0_6 = simplify(T0_5 * T5_6)
T0_G = simplify(T0_6 * T6_G)

print "Building total transformation matrix..."
# Correction for difference in gripper link orientation in URDF versus DH convention
# Transform for Body-fixed 180 deg rotation about z-axis
Tz = Matrix([[ cos(pi), -sin(pi), 0, 0],
             [ sin(pi),  cos(pi), 0, 0],
             [       0,        0, 1, 0],
             [       0,        0, 0, 1]])
# Transform for body-fixed -90 deg rotation about y-axis
Ty = Matrix([[  cos(-pi/2), 0, sin(-pi/2), 0],
             [           0, 1,          0, 0],
             [ -sin(-pi/2), 0, cos(-pi/2), 0],
             [           0, 0,          0, 1]])
# Apply the transforms to get a correction transformation matrix
T_corr = simplify(Tz * Ty)

# Total transformation matrix after orientation correction is applied
T_total = simplify(T0_G * T_corr)
#print('T_total: ', T_total)
print "Initialization COMPLETE!"

def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            #print(req.poses[x].position)
            #print(req.poses[x].orientation)
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Joint angle symbols
            theta1, theta2, theta3, theta4, theta5, theta6 = symbols('theta1:7')

            # Extract end-effector position and orientation from request
	        # px,py,pz = end-effector position
	        # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])
     
            # Create yaw, pitch and roll angle rotation matrices
            Ry = Matrix([[ cos(yaw), -sin(yaw), 0],
                         [ sin(yaw),  cos(yaw), 0],
                         [        0,         0, 1]])
            
            Rp = Matrix([[  cos(pitch), 0, sin(pitch)],
                         [           0, 1,          0],
                         [ -sin(pitch), 0, cos(pitch)]])
            
            Rr = Matrix([[ 1,         0,          0],
                         [ 0, cos(roll), -sin(roll)],
                         [ 0, sin(roll),  cos(roll)]])
                       
            #print(Rr, Rp, Ry)
            
            # Compute the rotation matrix Rrpy of the rotation between base_link and gripper_link
            Rrpy = Ry * Rp * Rr

            # Undo the orientation correction to the Rrpy to get the DH version
            Rrpy = simplify(Rrpy * T_corr[0:3,0:3].inv())
            #print(Rrpy)
            
            # Calculate joint angles using Geometric IK method
            
            # First, calculate the wrist center position using inverse position problem
            # Gripper length
            gl = (d6+d7).subs(s)
            
            # Wrist center position
            Wx = simplify(px - gl*Rrpy[0,2])
            Wy = simplify(py - gl*Rrpy[1,2])
            Wz = simplify(pz - gl*Rrpy[2,2])
            
            # Set r to XY distance of wrist center
            r = sqrt(Wx**2 + Wy**2)
            #print("W: ", Wx, Wy, Wz, " r: ", r)
            
            # Calculate theta1
            theta1 = atan2(Wy, Wx)
            #print("theta1: ", degrees(theta1))

            # Calculate theta2 using sympy to "solve" eq: a*cos(x) + b*sin(x) - c = 0
            eq = (d1-Wz)*cos(theta2) + (a1-r)*sin(theta2) - ((d4**2-(d1-Wz)**2-(a1-r)**2-(a2-a3)**2)/(2*(a2-a3)))
            eq = simplify(eq.subs(s))
            sol = solve(eq)
            #print("eq: ", eq, "   sol: ", sol)
            # Skip, if a sol does not exist
            if sol is None or len(sol)<1:
                continue
            # Otherwise, use the first valid value as theta2
            theta2 = sol[0]
            #print("theta2: ", degrees(theta2), degrees(sol[0]), degrees(sol[1]), eq)
            
            # Calculate theta3 by substituting theta2
            theta3 = atan2(((d1-Wz)+(a2-a3)*cos(theta2))/d4, -((a1-r)+(a2-a3)*sin(theta2))/d4) - theta2
            theta3 = simplify(theta3.subs(s))
            #print("theta3: ", degrees(theta3))
			
            # Inverse Orientation problem for finding theta 4,5 and 6
            
            # Substitute theta 1, 2 and 3 in T0_3 to get R0_3
            R0_3 = T0_3[0:3,0:3].evalf(subs={q1:theta1, q2:theta2, q3:theta3})

            # Determine R3_6 = inv(R0_3) * Rrpy, whose elements are functions of only theta 4, 5 and 6
            R3_6 = R0_3.inv() * Rrpy
            
            # Comparing individual elements to calculate theta 4,5 and 6
            
            # Calculate theta5 by solving cos(theta5) = R3_6[1,2]
            theta5 = simplify(acos(R3_6[1,2]))
            #print("theta4: ", degrees(theta4))
            
            # Calculate theta4 by solving sin(theta4)sin(theta5) = R3_6[2,2]
            theta4 = simplify(asin(R3_6[2,2]/sin(theta5)))
            #print("theta5: ", degrees(theta5))

            # Calculate theta6 by solving -sin(theta5)sin(theta6) = R3_6[1,1]
            theta6 = simplify(asin(R3_6[1,1]/-sin(theta5)))
            #print("theta6: ", degrees(theta6))
            
            # Populate response for the IK request theta1,theta2...,theta6 joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)
            print(len(joint_trajectory_list), len(req.poses), 'thetas: ', joint_trajectory_point.positions)
            
            # Calculate errors using forward kinematics, 
            # Substitute joint angles theta1-6 in T_total
            # Translation (last column) will be the end-effector position relative to base which is at origin
            T = T_total.evalf(subs={q1:theta1, q2:theta2, q3:theta3, q4:theta4, q5:theta5, q6:theta6})
            ex = px - T[0,3]
            ey = py - T[1,3]
            ez = pz - T[2,3]
            ecum = sqrt(ex**2 + ey**2 + ez**2)
            print ('errors: ', ex, ey, ez, 'cum: ', ecum)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK requests!"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
