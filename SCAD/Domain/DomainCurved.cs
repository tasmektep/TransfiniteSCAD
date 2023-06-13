using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SCAD
{
    class DomainCurved : Domain
    {

        ~DomainCurved() { }

        public override bool update()
        {
            int m = curves_.Count();
            List<Point2d> temp_vertices = new List<Point2d>();
            vertices_.resize(m);

            for (int i = 0; i < m; ++i)
            {
                var points = new Point3d[0];
                curves_[i].DivideByCount(resolution, true, out points);
                vertices_[i] = new Point2d(curves_[i % curves_.Count].PointAtStart.X, curves_[i % curves_.Count].PointAtStart.Y);
            }

            return base.update();
        }

        //protected override List<Point2d> ParametersImpl(int resolution)
        //{
        //    List<Point2d> parameters = new List<Point2d>(MeshSize(resolution));


        //    parameters.Add(center_);
        //    for (int j = 1; j <= resolution; ++j)
        //    {
        //        double u = (double)j / (double)resolution;
        //        for (int k = 0; k < n_; ++k)
        //            for (int i = 0; i < j; ++i)
        //            {
        //                double v = (double)i / (double)j;
        //                Point2d ep = vertices_.Prev(k) * (1.0 - v) + vertices_[k] * v;
        //                Point2d p = center_ * (1.0 - u) + ep * u;
        //                parameters.Add(p);
        //            }
        //    }

        //    return parameters;
        //}


        protected override void ComputeCenter() => center_ = new Point2d(0.0, 0.0);

    }
}