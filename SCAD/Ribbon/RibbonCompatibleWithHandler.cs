using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Media3D;

namespace SCAD
{
    class RibbonCompatibleWithHandler : RibbonCompatible
    {
        protected Vector3d central_;

        ~RibbonCompatibleWithHandler() { }

        public override void update()
        {
          base.update();
            Vector3d n = normal(0.5);
            if (!handler_initialized_)
            {
                handler_ = prev_tangent_ / prev_tangent_.Length + next_tangent_ / next_tangent_.Length;
                handler_ = handler_ - n * (handler_ * n);
                handler_.Unitize();
            }
            central_ = handler_ * (prev_tangent_.Length + next_tangent_.Length) / 2.0 * multiplier_;
        }
        public override Vector3d crossDerivative(double s)
        {
            Vector3d n = normal(s);
            Vector3d pt = prev_tangent_ - n * (prev_tangent_ * n);
            Vector3d ch = central_ - n * (central_ * n);
            Vector3d nt = next_tangent_ - n * (next_tangent_ * n);
            return pt * 2.0 * (s - 1.0) * (s - 0.5)
                + ch * -4.0 * s * (s - 1.0)
                + nt * 2.0 * s * (s - 0.5);
        }
    }
}
