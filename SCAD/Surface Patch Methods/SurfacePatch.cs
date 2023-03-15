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
    public class SurfacePatch
    {
        Domain dm;
        Parametrization prm;
        Blending_Method blending_method;

        /// Isocurve Mean value


        public SurfacePatch(Domain dm, Parametrization_Method method , Blending_Method blending_method = Blending_Method.Special_Side_Blending )
        {
            this.dm = dm;
            this.blending_method = blending_method; 
            if (method == Parametrization_Method.RadialDistanceFunction)
                prm = new Parametrization(method, dm);
            else if (method == Parametrization_Method.Harmonic_Pt)
                prm = new Parametrization(Parametrization_Method.Harmonic_Pt, dm);
        }

        /// <summary>
        /// kato
        /// curve üzerinde ribbonlara ihtiyaç var, eğer curve içinde gömülü değilse
        /// </summary>
        /// <param name="u"> u value </param>
        /// <param name="v"> v value </param>
        /// <param name="curves"> boundary curves </param>
        /// <returns></returns>
        /// 
        public Point3d Kato_Suv(double u, double v)
        {
            List<Curve> curves = dm.Curves;

            List<double> si = new List<double>();
            List<double> di = new List<double>();
            (si, di) = prm.GetPoint(u, v);

            BlendingFunctions blending = new BlendingFunctions(blending_method);
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
