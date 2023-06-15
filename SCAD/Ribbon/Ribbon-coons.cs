using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;
using static SCAD.Extensions;
namespace SCAD
{
    class RibbonCoons : Ribbon
    {
        protected NurbsCurve left_, right_, top_;
        protected Point3d bl_, br_, tl_, tr_;
        ~RibbonCoons() { }
        public override void Update()
        {
            left_ = prev_.Curve();
            right_ = next_.Curve();

            var left2 = prev_.prev_.Curve();
            var right2 = prev_.prev_.Curve();
            //var left2 = dynamic_cast<RibbonCoons*>(prev_.lock ().get())->prev_.lock ()->curve();
            //var right2 = dynamic_cast<RibbonCoons*>(next_.lock ().get())->next_.lock ()->curve();
            var p1 = left_.PointAt(0.0);
            var q1 = right_.PointAt(1.0);
            if (right2 == left_)          // 3-sided
                top_ = NurbsCurve.Create(false, 0, new List<Point3d>() { p1 }); // 0-degree Bezier curve

            else
            {
                Vector3d[] der;
                der = left2.DerivativeAt(1.0, 1);
                var p2 = p1 - der[1] / 3.0;
                der = right2.DerivativeAt(0.0, 1);
                var q2 = q1 + der[1] / 3.0;
                var lp  = new List<Point3d>() { q1, q2, p2, p1 };
                top_ = NurbsCurve.Create(false, 1, lp); // 0-degree Bezier curve
            }
            bl_ = curve_.PointAt(0.0);
            br_ = curve_.PointAt(1.0);
            tl_ = left_.PointAt(0.0);
            tr_ = right_.PointAt(1.0);
            base.Update();
        }
        public override Vector3d CrossDerivative(double s)
        {
            // TODO
            return new Vector3d();
        }
        public override Point3d Eval(Point2d sd)
        {
            double s = inrange(0, sd[0], 1), d = inrange(0, sd[1], 1),
             s1 = inrange(0, 1 - s, 1), d1 = inrange(0, 1 - d, 1);
            var p1 = curve_.PointAt(s) * d1 + top_.PointAt(s1) * d;
            var p2 = left_.PointAt(d1) * s1 + right_.PointAt(d) * s;
            var p12 = (bl_ * s1 + br_ * s) * d1 + (tl_ * s1 + tr_ * s) * d;
            var vec = p1 + p2 - p12;
            return new Point3d(vec.X, vec.Y, vec.Z);
        }


    }
}
