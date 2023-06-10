using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SCAD
{
    internal class DomainConcave: Domain
    {
        ~DomainConcave() { }

        public override bool update()
        {
            int m = curves_.Count();


            vertices_.resize(m);
            for (int i = 0; i < m; ++i)
            {

                vertices_[i] = new Point2d(curves_[i % curves_.Count].PointAtStart.X, curves_[i % curves_.Count].PointAtStart.Y);
            }


            return base.update();
        }
        protected override void ComputeCenter() => center_ = new Point2d(0.0, 0.0);


    }
}
