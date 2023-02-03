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
        private Point3d Kato_Suv(double u, double v, List<Curve> curves)
        {
            new Domain(curves, out List<Curve> domaincurve);
            Parametrization prm = new Parametrization(Parametrization_Method.Harmonic,domaincurve);
            List<double> si = new List<double>();
            List<double> di = new List<double>();
            (si, di) = prm.GetPoint(u, v);

            List<double> d_i = new List<double>();
            for (int i = 0; i < si.Count; i++)
            {
                d_i.Add(di[i]);
            }
            BlendingFunctions blending = new BlendingFunctions(Blending_Method.Special_Side_Blending);
            List<double> Value = blending.GetBlending(d_i);

            Point3d r_sum = new Point3d();
            for (int i = 0; i < curves.Count; i++)
            {
                double s = curves[i].Domain.Min + si[i] * (curves[i].Domain.Max - curves[i].Domain.Min);

                Vector3d T = new Vector3d(); // Ribbon Vector EKSİKKKKK
                Point3d r = curves[i].PointAt(s) + (di[i] * T);

                r_sum += r * Value[i];
            }

            return r_sum;
        }

    }
}
