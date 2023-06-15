using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SCAD
{
    class DomainRegular : Domain
    {
        ~DomainRegular() { }

        public override bool Update()
        {
            int m = curves_.Count();
            if (n_ == m)
                return false;

            if (m == 4)
            {
                vertices_ = new List<Point2d>() { new Point2d(1, 1), new Point2d(-1, 1), new Point2d(-1, -1), new Point2d(1, -1) };
                return base.Update();
            }

            double alpha = 2.0 * M_PI / m;
            vertices_.resize(m);
            for (int i = 0; i < m; ++i)
            {
                int j = (i + 1) % curves_.Count;
                vertices_[i] = new Point2d(Math.Cos(alpha * i), Math.Sin(alpha * i));

            }
            //vertices_[i] = Point2D(std::cos(alpha * i), std::sin(alpha * i));


            return base.Update();
        }
        protected override void ComputeCenter() => center_ = new Point2d(0.0, 0.0);
    }
}