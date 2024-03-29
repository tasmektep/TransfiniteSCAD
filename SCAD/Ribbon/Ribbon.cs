﻿using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
using System.Security.Cryptography;
using System.Windows.Media.Media3D;

namespace SCAD
{
    public class Ribbon
    {
        protected NurbsCurve curve_;
        public Ribbon prev_, next_;
        protected RMF rmf_ = new RMF();
        protected Vector3d normal_fence_ = new Vector3d();
        protected Vector3d handler_;
        protected double multiplier_;
        protected bool handler_initialized_;

        public NurbsCurve curve()
        {
            return curve_;
        }

        public void setCurve(NurbsCurve curve)
        {
            curve_ = curve;
        }

        public void setNeighbors(Ribbon prev, Ribbon next)
        {
            prev_ = prev;
            next_ = next;
        }

        public double multiplier()
        {
            return multiplier_;
        }

        public void setMultiplier(double m)
        {
            multiplier_ = m;
        }

        public Vector3d handler()
        {
            if (handler_initialized_)
                return handler_;
            return new Vector3d();
        }

        public void setHandler(Vector3d h)
        {
            handler_ = h;
            handler_.Unitize();
            handler_initialized_ = true;
        }


        public void overrideNormalFence(Vector3d fence)
        {
            normal_fence_ = fence;
        }
        public void reset()
        {
            multiplier_ = 1.0;
            handler_initialized_ = false;
        }
        public virtual void update()
        {
            Vector3d normal;
            Vector3d[] der;
            rmf_.setCurve(curve_);
            der = prev_.curve_.DerivativeAt(1.0, 1);
            normal = der[1];
            der = curve_.DerivativeAt(0.0, 1);
            normal = Vector3d.CrossProduct(normal, der[1]);
            normal.Unitize();
            rmf_.setStart(normal);

            der = curve_.DerivativeAt(1.0, 1);
            normal = der[1];
            der = next_.curve_.DerivativeAt(0.0, 1);
            normal = Vector3d.CrossProduct(normal, der[1]);
            normal.Unitize();
            rmf_.setEnd(normal);

            rmf_.update();
        }

        public virtual Vector3d crossDerivative(double s)
        { return new Vector3d(); }
        public virtual Point3d eval(Point2d sd)
        {
            return curve_.PointAt(sd[0]) + crossDerivative(sd[0]) * sd[1];
        }
        public Vector3d normal(double s)
        {
            //if (normal_fence_)
            //    return normal_fence_->operator()(s);
            return rmf_.eval(s);
        }
    }
}
