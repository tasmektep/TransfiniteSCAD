using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper;
using Rhino.Geometry;
using static SCAD.Extensions;

namespace SCAD
{
    class SurfacePatch
    {

        ////
        ///Blending function bana di değeri girdisi verildiğinde n (yani curve sayısı) kadar value değeri vermeli.

        /// <summary>
        /// kato
        /// curve üzerinde ribbonlara ihtiyaç var, eğer curve içinde gömülü değilse
        /// </summary>
        /// <param name="u"> u value </param>
        /// <param name="v"> v value </param>
        /// <param name="curves"> boundary curves </param>
        /// <returns></returns>
        public Point3d Kato_Suv(double u, double v, List<Curve> curves)
        {
            var dm = new domain();
            dm.setSides(curves);
            dm.Update();

            var mesh = dm.MeshTopology(30);
            var uvs = dm.parameters(30);
            var vs = dm.vertices_;


            List<Line> lines = new List<Line>();
            for (int i = 0; i < 6; i++)
            {
                lines.Add(new Line(new Point3d(vs[i].X, vs[i].Y, 0), new Point3d(vs[(i + 1) % 6].X, vs[(i + 1) % 6].Y, 0)));
            }
            //new Domain(curves, out List<Curve> domaincurve);
            //new Domain(curves, out List<Curve> domaincurve);
            //Parametrization prm = new Parametrization(Parametrization_Method.Harmonic,domaincurve);
            //List<double> si = new List<double>();
            //List<double> di = new List<double>();
            //(si, di) = prm.GetPoint(u, v);

            Parametrization prm = new Parametrization(Parametrization_Method.RadialDistanceFunction, lines);
            List<double> si = new List<double>();
            List<double> di = new List<double>();
            (si, di) = prm.GetPoint(u, v);

            BlendingFunctions blending = new BlendingFunctions(Blending_Method.Special_Side_Blending);
            List<double> Value = blending.GetBlending(di);

            Point3d r_sum = new Point3d();
            for (int i = 0; i < curves.Count; i++)
            {
                double s = curves[i].Domain.Min + si[i] * (curves[i].Domain.Max - curves[i].Domain.Min);
                Vector3d crossproduct = Vector3d.CrossProduct(curves[i].TangentAt(s), curves[i].CurvatureAt(s));

                Vector3d T = crossproduct; // Ribbon vector
                Point3d r = curves[i].PointAt(s) + (di[i] * T);

                r_sum += r * Value[i];
            }

            return r_sum;
        }

    }
}
