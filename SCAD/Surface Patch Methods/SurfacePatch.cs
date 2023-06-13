using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper;
using Rhino.Geometry;
using Rhino.Input.Custom;
using static SCAD.Extensions;

namespace SCAD
{
    public class SurfacePatch
    {
        Domain dm;
        Parametrization prm;
        Blending_Method blending_method;
        public List<Point3d> centers = new List<Point3d>();
        public List<Vector3d> vectors = new List<Vector3d>();
        public List<Plane> planes = new List<Plane>();
        public List<Ribbon> ribbons = new List<Ribbon>();
        /// Isocurve Mean value

        public SurfacePatch()
        { }

        public Parametrization GetParametrization { get {return prm; } }
        public SurfacePatch(Domain dm, List<Ribbon> ribbons, Parametrization_Method method, Blending_Method blending_method = Blending_Method.Special_Side_Blending)
        {
            this.dm = dm;
            this.blending_method = blending_method;
            if (method == Parametrization_Method.RadialDistanceFunction)
                prm = new Parametrization(method, dm);
            else if (method == Parametrization_Method.Harmonic_Pt)
                prm = new Parametrization(Parametrization_Method.Harmonic_Pt, dm);
            this.ribbons = ribbons;
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
                Plane VecPlane = new Plane();
                //double s = curves[i].Domain.Min + si[i] * (curves[i].Domain.Max - curves[i].Domain.Min);
                curves[i].Domain = new Interval(0.0, 1.0);
                Vector3d crossproduct = Vector3d.CrossProduct(curves[i].TangentAt(si[i]), curves[i].CurvatureAt(si[i]));
                //Vector3d T = (ribbons[i].eval(new Point2d(si[i], di[i]))- curves[i].PointAt(si[i])); // Ribbon vector      
                Vector3d T = (ribbons[i].eval(new Point2d(si[i], 1.0))- curves[i].PointAt(si[i])); // Ribbon vector      
                //T.Unitize();
                planes.Add(VecPlane);
                Point3d r = curves[i].PointAt(si[i]) + (di[i] * T);
                vectors.Add(ribbons[i].eval(new Point2d(si[i], di[i])) - curves[i].PointAt(si[i]));
                centers.Add(curves[i].PointAt(si[i]));
                r_sum += r * Value[i];
            }

            return r_sum;
        }
    }
}
