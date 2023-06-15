using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SCAD
{
    internal class DomainCurved : Domain
    {
        ~DomainCurved() { }

        public override bool Update()
        {
            int m = curves_.Count();
            var temp_vertices = new List<Point3d>();
            var temp_vertices2d = new List<Point2d>();

            int div_count =3;
            vertices_.resize(m);
            for (int i = 0; i < m; ++i)
            {
                //if (i == 3)
                //{
                //    curves_[i].DivideByCount(1, true, out Point3d[] tempts);
                //    temp_vertices.AddRange(tempts);
                //}
                //else
                //{
                    curves_[i].DivideByCount(div_count, true, out Point3d[] tempts);
                    temp_vertices.AddRange(tempts);
                //}
            }

            foreach (var temp_vertice in temp_vertices)
            {
                temp_vertices2d.Add(new Point2d(temp_vertice.X, temp_vertice.Y));
            }
            vertices_ = temp_vertices2d;
            return base.Update();
        }
        protected override void ComputeCenter() => center_ = new Point2d(0.0, 0.0);


    }
}
