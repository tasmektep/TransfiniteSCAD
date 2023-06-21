using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Net.Mime.MediaTypeNames;

namespace SCAD
{
    internal class RibbonScad : Ribbon
    {
        protected Vector3d prev_tangent_, next_tangent_;
        ~RibbonScad() { }
        public override void Update()
        {
            Vector3d[] der;
            der = prev_.Curve().DerivativeAt(1.0, 1);
            prev_tangent_ = -der[1];
            der = next_.Curve().DerivativeAt(0.0, 1);
            next_tangent_ = der[1];
            base.Update();
        }

        public Vector3d CrossDerivative(double s, Vector3d sphereNormal)
        {
            Vector3d n = sphereNormal;
            Vector3d pt = prev_tangent_ - n * (prev_tangent_ * n);
            Vector3d nt = next_tangent_ - n * (next_tangent_ * n);
            return pt * (1.0 - s) + nt * s;
        }

    }
}
