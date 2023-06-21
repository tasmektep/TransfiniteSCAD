using Rhino.Geometry;

namespace SCAD
{
    public enum Ribbon_Method
    {
        //
        // Summary:
        //      Uses Side Blending functions
        Ribbon_Compatible = 0,
        //
        // Summary:
        //      Uses Corner Blending functions
        Ribbon_CompatiablewithHandler,
        //
        // Summary:
        //      Uses Special Side Blending functions
        Ribbon_Coons
    }

    public class Ribbon
    {

        protected NurbsCurve curve_;
        public Ribbon prev_, next_;
        protected RMF rmf_ = new RMF();
        protected Vector3d normal_fence_ = new Vector3d();
        protected Vector3d handler_;
        protected double multiplier_;
        protected bool handler_initialized_;

        public NurbsCurve Curve()
        {
            return curve_;
        }

        public void SetCurve(NurbsCurve curve)
        {
            curve_ = curve;
        }

        public void SetNeighbors(Ribbon prev, Ribbon next)
        {
            prev_ = prev;
            next_ = next;
        }

        public double Multiplier()
        {
            return multiplier_;
        }

        public void SetMultiplier(double m)
        {
            multiplier_ = m;
        }

        public Vector3d Handler()
        {
            if (handler_initialized_)
                return handler_;
            return new Vector3d();
        }

        public void SetHandler(Vector3d h)
        {
            handler_ = h;
            handler_.Unitize();
            handler_initialized_ = true;
        }


        public void OverrideNormalFence(Vector3d fence)
        {
            normal_fence_ = fence;
        }
        public void Reset()
        {
            multiplier_ = 1.0;
            handler_initialized_ = false;
        }
        public virtual void Update()
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

        public virtual Vector3d CrossDerivative(double s)
        { 
            return new Vector3d(); 
        }

        public virtual Vector3d CrossDerivative(double s, Vector3d norm)
        {
            return new Vector3d();
        }

        public virtual Point3d Eval(Point2d sd)
        {
            return curve_.PointAt(sd[0]) + CrossDerivative(sd[0]) * sd[1];
        }
        public Vector3d Normal(double s)
        {
            //if (normal_fence_)
            //    return normal_fence_->operator()(s);
            return rmf_.eval(s);
        }
    }
}
