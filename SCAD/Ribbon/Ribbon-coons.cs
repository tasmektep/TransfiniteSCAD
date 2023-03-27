using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;

namespace SCAD
{
    class RibbonCoons : Ribbon
    {
        protected NurbsCurve left_, right_, top_;
        protected Point3d bl_, br_, tl_, tr_;
        ~RibbonCoons() { }
        public override void update()
        {
            left_ = prev_.curve();
            right_ = next_.curve();

            var left2 = prev_.prev_.curve();
            var right2 = prev_.prev_.curve();
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
                top_ = NurbsCurve.Create(false, 0, new List<Point3d>() { q1, q2, p2, p1 }); // 0-degree Bezier curve
            }
            bl_ = curve_.PointAt(0.0);
            br_ = curve_.PointAt(1.0);
            tl_ = left_.PointAt(0.0);
            tr_ = right_.PointAt(1.0);
            base.update();
        }
        //      public virtual Vector3D crossDerivative(double s) const override;
        //public virtual Point3D eval(const Point2D &sd) const override;


    }
}
