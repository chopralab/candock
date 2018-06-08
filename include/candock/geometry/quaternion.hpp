#ifndef QUATERNION_H
#define QUATERNION_H
#include "candock/helper/debug.hpp"
#include "candock/geometry/geometry.hpp"
#include "candock/geometry/matrix.hpp"
#include <tuple>
#include <functional>

#include <cmath>

/*
 * Taken from : http://www.stanford.edu/~acoates/quaternion.h
 * 
 */ 
namespace geometry {
	class Quaternion {
		double mData[4];
	public:
		Quaternion() {
			mData[0] = mData[1] = mData[2] = 0;
			mData[3] = 1;
		}
		Quaternion(const Vector3& v, double w) {
			//~ mData[0] = v.element(0,0);
			//~ mData[1] = v.element(1,0);
			//~ mData[2] = v.element(2,0);
			mData[0] = v.x();
			mData[1] = v.y();
			mData[2] = v.z();
			mData[3] = w;
		}
		Quaternion(double x, double y, double z, double w) {
			mData[0] = x;
			mData[1] = y;
			mData[2] = z;
			mData[3] = w;
		}		
		double x() const { return mData[0]; }
		double y() const { return mData[1]; }
		double z() const { return mData[2]; }
		double w() const { return real(); }
		
		Vector3 complex() const { return Vector3(mData); }
		//~ void complex(const Vector3& c) { mData[0] = c[0]; mData[1] = c[1];  mData[2] = c[2]; }
		void complex(const Vector3& c) { mData[0] = c.x(); mData[1] = c.y();  mData[2] = c.z(); }
		
		double real() const { return mData[3]; }
		void real(double r) { mData[3] = r; }

		Quaternion conjugate(void) const {
			return Quaternion(-complex(), real());
		}
		
		/** 
		* @brief Computes the inverse of this quaternion.
		*
		* @note This is a general inverse.  If you know a priori
		* that you're using a unit quaternion (i.e., norm() == 1),
		* it will be significantly faster to use conjugate() instead.
		* 
		* @return The quaternion q such that q * (*this) == (*this) * q
		* == [ 0 0 0 1 ]<sup>T</sup>.
		*/
		Quaternion inverse(void) const {
			return conjugate() / norm();
		}

		/** 
		* @brief Computes the product of this quaternion with the
		* quaternion 'rhs'.
		*
		* @param rhs The right-hand-side of the product operation.
		*
		* @return The quaternion product (*this) x @p rhs.
		*/
		Quaternion product(const Quaternion& rhs) const {
			return Quaternion(y()*rhs.z() - z()*rhs.y() + x()*rhs.w() + w()*rhs.x(),
							z()*rhs.x() - x()*rhs.z() + y()*rhs.w() + w()*rhs.y(),
							x()*rhs.y() - y()*rhs.x() + z()*rhs.w() + w()*rhs.z(),
							w()*rhs.w() - x()*rhs.x() - y()*rhs.y() - z()*rhs.z());
		}
		
		/**
		* @brief Quaternion product operator.
		*
		* The result is a quaternion such that:
		*
		* result.real() = (*this).real() * rhs.real() -
		* (*this).complex().dot(rhs.complex());
		*
		* and:
		*
		* result.complex() = rhs.complex() * (*this).real
		* + (*this).complex() * rhs.real()
		* - (*this).complex().cross(rhs.complex());
		*
		* @return The quaternion product (*this) x rhs.
		*/
		Quaternion operator*(const Quaternion& rhs) const {
			return product(rhs);
		}
		/**
		* @brief Quaternion scalar product operator.
		* @param s A scalar by which to multiply all components
		* of this quaternion.
		* @return The quaternion (*this) * s.
		*/
		Quaternion operator*(double s) const {
			return Quaternion(complex()*s, real()*s);
		}
		/**
		* @brief Produces the sum of this quaternion and rhs.
		*/
		Quaternion operator+(const Quaternion& rhs) const {
			return Quaternion(x()+rhs.x(), y()+rhs.y(), z()+rhs.z(), w()+rhs.w());
		}
		
		/**
		* @brief Produces the difference of this quaternion and rhs.
		*/
		Quaternion operator-(const Quaternion& rhs) const {
			return Quaternion(x()-rhs.x(), y()-rhs.y(), z()-rhs.z(), w()-rhs.w());
		}
		
		/**
		* @brief Unary negation.
		*/
		Quaternion operator-() const {
			return Quaternion(-x(), -y(), -z(), -w());
		}
		/**
		* @brief Quaternion scalar division operator.
		* @param s A scalar by which to divide all components
		* of this quaternion.
		* @return The quaternion (*this) / s.
		*/
		Quaternion operator/(double s) const {
			//~ if (s == 0) std::clog << "Dividing quaternion by 0." << std::endl;
			if (s == 0) throw Error("die : Dividing quaternion by 0.");
			return Quaternion(complex()/s, real()/s);
		}
		/**
		* @brief Returns the norm ("magnitude") of the quaternion.
		* @return The 2-norm of [ w(), x(), y(), z() ]<sup>T</sup>.
		*/
		double norm() const { return sqrt(mData[0]*mData[0]+mData[1]*mData[1]+
							  mData[2]*mData[2]+mData[3]*mData[3]); }
		
		/**
		* @brief Computes the rotation matrix represented by a unit
		* quaternion.
		*
		* @note This does not check that this quaternion is normalized.
		* It formulaically returns the matrix, which will not be a
		* rotation if the quaternion is non-unit.
		*/
		Matrix rotationMatrix() const {
			double m[9] = {
				1-2*y()*y()-2*z()*z(), 2*x()*y() - 2*z()*w(), 2*x()*z() + 2*y()*w(),
				2*x()*y() + 2*z()*w(), 1-2*x()*x()-2*z()*z(), 2*y()*z() - 2*x()*w(),
				2*x()*z() - 2*y()*w(), 2*y()*z() + 2*x()*w(), 1-2*x()*x()-2*y()*y()
			};
			return Matrix(m);
		}
		/**
		* @brief Returns a vector rotated by this quaternion.
		*
		* Functionally equivalent to:  (rotationMatrix() * v)
		* or (q * Quaternion(0, v) * q.inverse()).
		*
		* @warning conjugate() is used instead of inverse() for better
		* performance, when this quaternion must be normalized.
		*/
		Vector3 rotatedVector(const Vector3& v) const {
			return (((*this) * Quaternion(v, 0)) * conjugate()).complex();
		}
		/**
		* @brief Computes the quaternion that is equivalent to a given
		* euler angle rotation.
		* @param euler A 3-vector in order:  roll-pitch-yaw.
		*/
		void euler(const Vector3& euler) {
			//~ double c1 = cos(euler[2] * 0.5);
			//~ double c2 = cos(euler[1] * 0.5);
			//~ double c3 = cos(euler[0] * 0.5);
			//~ double s1 = sin(euler[2] * 0.5);
			//~ double s2 = sin(euler[1] * 0.5);
			//~ double s3 = sin(euler[0] * 0.5);
			double c1 = cos(euler.z() * 0.5);
			double c2 = cos(euler.y() * 0.5);
			double c3 = cos(euler.x() * 0.5);
			double s1 = sin(euler.z() * 0.5);
			double s2 = sin(euler.y() * 0.5);
			double s3 = sin(euler.x() * 0.5);
			
			mData[0] = c1*c2*s3 - s1*s2*c3;
			mData[1] = c1*s2*c3 + s1*c2*s3;
			mData[2] = s1*c2*c3 - c1*s2*s3;
			mData[3] = c1*c2*c3 + s1*s2*s3;
		}
		
		/** @brief Returns an equivalent euler angle representation of
		* this quaternion.
		* @return Euler angles in roll-pitch-yaw order.
		*/
		Vector3 euler(void) const {
			Vector3 euler;
			const static double PI_OVER_2 = M_PI_2;
			const static double EPSILON = 1e-10;
			double sqw, sqx, sqy, sqz;
			
			// quick conversion to Euler angles to give tilt to user
			sqw = mData[3]*mData[3];
			sqx = mData[0]*mData[0];
			sqy = mData[1]*mData[1];
			sqz = mData[2]*mData[2];
			
			//~ euler[1] = asin(2.0 * (mData[3]*mData[1] - mData[0]*mData[2]));
			euler.set_y(asin(2.0 * (mData[3]*mData[1] - mData[0]*mData[2])));
			//~ if (PI_OVER_2 - fabs(euler[1]) > EPSILON) {
			if (PI_OVER_2 - fabs(euler.y()) > EPSILON) {
				//~ euler[2] = atan2(2.0 * (mData[0]*mData[1] + mData[3]*mData[2]),
				euler.set_z(atan2(2.0 * (mData[0]*mData[1] + mData[3]*mData[2]),
				sqx - sqy - sqz + sqw));
				//~ euler[0] = atan2(2.0 * (mData[3]*mData[0] + mData[1]*mData[2]),
				euler.set_x(atan2(2.0 * (mData[3]*mData[0] + mData[1]*mData[2]),
				sqw - sqx - sqy + sqz));
			} else {
				// compute heading from local 'down' vector
				//~ euler[2] = atan2(2*mData[1]*mData[2] - 2*mData[0]*mData[3],
				//~ 2*mData[0]*mData[2] + 2*mData[1]*mData[3]);
				euler.set_z(atan2(2*mData[1]*mData[2] - 2*mData[0]*mData[3],
				2*mData[0]*mData[2] + 2*mData[1]*mData[3]));
				//~ euler[0] = 0.0;
				euler.set_x(0.0);
				
				// If facing down, reverse yaw
				//~ if (euler[1] < 0)
				if (euler.y() < 0)
				//~ euler[2] = M_PI - euler[2];
				euler.set_z(M_PI - euler.z());
			}
			return euler;
		}
		
	};
}

#endif
