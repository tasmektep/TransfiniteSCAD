using System;
using System.Collections.Generic;
using System.Diagnostics.PerformanceData;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Markup;
using System.Windows.Media.Media3D;
using Grasshopper;
using Grasshopper.GUI;
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
        private int n_;
        const double epsilon = 1.0e-8;
        bool use_gamma_ = true;
        private List<CornerData> corner_data_ = new List<CornerData>();
        List<NurbsCurve> hermits = new List<NurbsCurve>();
        List<Point3d> point3Ds = new List<Point3d>();
        List<Point2d> uvss = new List<Point2d>();
        double eps = 0.0000001;

        /// Isocurve Mean value

        public SurfacePatch()
        { }

        public List<NurbsCurve> GetHermits { get { return hermits; } }

        public List<Point3d> GetPoints { get { return point3Ds; } }

        public List<Point2d> GetUv { get { return uvss; } }


        public Parametrization GetParametrization { get { return prm; } }

        public SurfacePatch(Domain dm, List<Ribbon> ribbons, Parametrization_Method method, Blending_Method blending_method = Blending_Method.Special_Side_Blending)
        {
            this.dm = dm;
            this.blending_method = blending_method;
            if (method == Parametrization_Method.RadialDistanceFunction)
                prm = new Parametrization(method, dm);
            else if (method == Parametrization_Method.Harmonic_Pt)
                prm = new Parametrization(Parametrization_Method.Harmonic_Pt, dm);
            this.ribbons = ribbons;
            n_ = this.dm.Size;

        }

        public SurfacePatch(Domain dm, List<Ribbon> ribbons, List<CornerData> corner_data_, Parametrization_Method method, Blending_Method blending_method = Blending_Method.Special_Side_Blending)
        {
            this.corner_data_ = corner_data_;
            this.dm = dm;
            this.blending_method = blending_method;
            if (method == Parametrization_Method.RadialDistanceFunction)
                prm = new Parametrization(method, dm);
            else if (method == Parametrization_Method.Harmonic_Pt)
                prm = new Parametrization(Parametrization_Method.Harmonic_Pt, dm);
            this.ribbons = ribbons;
            n_ = this.dm.Size;
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

        public Point3d GetSurfacePoint(double u, double v)
        {
            //var pt = Kato_Suv(u, v);
            //var pt = GeneralizedCoon(u, v);
            var pt = Kato__hermit_Suv(u, v);

            return pt;
        }

        public Point3d Kato_Suv(double u, double v)
        {
            List<Curve> curves = dm.Curves;
            List<Curve> secondcurve = dm.TCurves;///***** (bunun içi ikinci 6 curve ile doldurulacak)
            Brep sphere = dm.Sphere;

            List<double> si, di;
            (si, di) = prm.GetPoint(u, v);

            BlendingFunctions blending = new BlendingFunctions(blending_method);
            List<double> Value = blending.GetBlending(di);

            Point3d r_sum = new Point3d();

            for (int i = 0; i < curves.Count; i++)
            {
                Plane VecPlane = new Plane();
                //double s = curves[i].Domain.Min + si[i] * (curves[i].Domain.Max - curves[i].Domain.Min);
                curves[i].Domain = new Interval(0.0, 1.0);
                //Vector3d crossproduct = Vector3d.CrossProduct(-curves[i].CurvatureAt(si[i]), curves[i].TangentAt(si[i]));

                //Vector3d tressholdproduct = secondcurve[i].PointAt(si[i]) - curves[i].PointAt(si[i]); ///*****
                //tressholdproduct.Unitize();///*****

                sphere.Surfaces[0].ClosestPoint(curves[i].PointAt(si[i]), out double longparam, out double latparam);
                Vector3d spherevec = Vector3d.CrossProduct(sphere.Surfaces[0].NormalAt(longparam, latparam), curves[i].TangentAt(si[i]));
                spherevec.Unitize();
                //spherevec = spherevec * 2.0;

                //Vector3d T = (ribbons[i].eval(new Point2d(si[i], di[i]))- curves[i].PointAt(si[i])); // Ribbon vector      
                Vector3d T = (ribbons[i].Eval(new Point2d(si[i], 1.0)) - curves[i].PointAt(si[i])); // Ribbon vector      
                T = new Vector3d(T.X, T.Y, 0.2);
                if (i == 0)
                    T = new Vector3d(0, 1, 0.2);
                if (i == 1)
                    T = new Vector3d(-1, 0, 0.2);
                if (i == 2)
                {
                    T = (1 - si[i]) * new Vector3d(0, -1, 0.2) + si[i] * new Vector3d(-1, 0, 0.2);
                }
                //T = new Vector3d(0, -1, 0.2);
                if (i == 3)
                    T = new Vector3d(-1, 0, 0.2);
                if (i == 4)
                    T = new Vector3d(0, -1, 0.2);
                if (i == 5)
                    T = new Vector3d(1, 0, 0.2);
                T.Unitize();
                planes.Add(VecPlane);

                //Point3d r = curves[i].PointAt(si[i]) + (di[i] * T);
                //Point3d r = curves[i].PointAt(si[i]) + (di[i] * tressholdproduct);//(di[i] * T) ///*****
                Point3d r = curves[i].PointAt(si[i]) + (di[i] * spherevec);//(di[i] * T)
                vectors.Add(ribbons[i].Eval(new Point2d(si[i], di[i])) - curves[i].PointAt(si[i]));
                centers.Add(curves[i].PointAt(si[i]));
                r_sum += r * Value[i];
            }

            return r_sum;
        }

        public Point3d Kato__hermit_Suv(double u, double v)
        {
            List<Curve> curves = dm.Curves;
            List<Curve> secondcurve = dm.TCurves;///***** (bunun içi ikinci 6 curve ile doldurulacak)
            Brep sphere = dm.Sphere;


            List<double> si, di;
            (si, di) = prm.GetPoint(u, v);

            BlendingFunctions blending = new BlendingFunctions(blending_method);

            List<double> Value = blending.GetBlending(di);

            Point3d r_sum = new Point3d();

            uvss.Add(new Point2d(u, v));

            for (int i = 0; i < curves.Count; i++)
            {
                var twin = (i + curves.Count / 2) % curves.Count;
                var s = si[i];
                double sTwin = 1 - s;

                if (s < eps && s>0)
                { s = 0; sTwin = 1; }
                if (s > 1)
                { s = 1; sTwin = 0; }

                var sphereVec1 = CrossProductVec(sphere, curves[i], s);
                var sphereVec2 = CrossProductVec(sphere, curves[twin], sTwin);
                var pts = new List<Point3d>() { curves[i].PointAt(s), curves[twin].PointAt(sTwin) };

                var hermit = NurbsCurve.CreateHSpline(pts, sphereVec1, -sphereVec2).ToNurbsCurve();
                hermit.Domain = new Interval(0.0, 1.0);
                Point3d r = hermit.PointAt(di[i]);

                point3Ds.Add(r);
                r_sum += r * Value[i];
                hermits.Add(hermit);

            }
            point3Ds.Add(r_sum);
            return r_sum;
        }

        public Vector3d CrossProductVec(Brep sphere_, Curve crv, double si)
        {
            sphere_.Surfaces[0].ClosestPoint(crv.PointAt(si), out double longparam1, out double latparam1);
            Vector3d sphereVec1 = Vector3d.CrossProduct(sphere_.Surfaces[0].NormalAt(longparam1, latparam1), crv.TangentAt(si));
            sphereVec1.Unitize();
            return sphereVec1;
        }


        public Point3d GeneralizedCoon(double u, double v)
        {
            List<Curve> curves = dm.Curves;

            List<double> si, di;
            (si, di) = prm.GetPoint(u, v);

            List<double> blends = GetBlending(si, di);
            Point3d p = new Point3d();
            for (int i = 0; i < curves.Count; ++i)
            {
                double s = si[i], d = di[i], s1 = si.Next(i);
                p += SideInterpolant(i, s, d) * (blends[i] + blends.Prev(i));
                p -= ((Vector3d)(CornerCorrection(i, 1.0 - s, s1) * blends[i]));
            }
            return p;
        }

        protected Point3d CornerCorrection(int i, double s1, double s2)
        {
            // Assumes that both s1 and s2 are 0 at the corner,
            // s1 increases towards corner (i-1), and s2 towards corner (i+1).
            s1 = inrange(0, Gamma(s1), 1);
            s2 = inrange(0, Gamma(s2), 1);
            return corner_data_[i].point
              + corner_data_[i].tangent1 * s1
              + corner_data_[i].tangent2 * s2
              + RationalTwist(s1, s2, corner_data_[i].twist2, corner_data_[i].twist1) * s1 * s2;
        }

        protected Point3d SideInterpolant(int i, double si, double di)
        {
            si = inrange(0, si, 1);
            di = Math.Max(Gamma(di), 0.0);
            return ribbons[i].Eval(new Point2d(si, di));
        }

        private List<double> GetBlending(List<double> si, List<double> di)
        {
            List<Point2d> sds = new List<Point2d>();
            for (int i = 0; i < di.Count; i++)
            {
                sds.Add(new Point2d(si[i], di[i]));
            }
            return BlendCorner(sds);
        }

        private double Gamma(double d)
        {
            if (use_gamma_)
                return d / (2.0 * d + 1.0);
            return d;
        }

        private static Vector3d RationalTwist(double u, double v, Vector3d f, Vector3d g)
        {
            if (Math.Abs(u + v) < epsilon)
                return new Vector3d(0, 0, 0);
            return (f * u + g * v) / (u + v);
        }

        protected List<double> BlendCorner(List<Point2d> sds)
        {
            List<double> blf = new List<double>();
            blf.Capacity = n_;

            int close_to_boundary = 0;
            foreach (var sd in sds)
            {
                if (sd[1] < epsilon)
                    ++close_to_boundary;
            }

            if (close_to_boundary > 0)
            {
                for (int i = 0; i < n_; ++i)
                {
                    int ip = Next(i);
                    if (close_to_boundary > 1)
                        blf.Add(sds[i][1] < epsilon && sds[ip][1] < epsilon ? 1.0 : 0.0);
                    else if (sds[i][1] < epsilon)
                    {
                        double tmp = Math.Pow(sds[ip][1], -2);
                        blf.Add(tmp / (tmp + Math.Pow(sds[Prev(i)][1], -2)));
                    }
                    else if (sds[ip][1] < epsilon)
                    {
                        double tmp = Math.Pow(sds[i][1], -2);
                        blf.Add(tmp / (tmp + Math.Pow(sds[Next(ip)][1], -2)));
                    }
                    else
                        blf.Add(0.0);
                }
            }
            else
            {
                double denominator = 0.0;
                for (int i = 0; i < n_; ++i)
                {
                    blf.Add(Math.Pow(sds[i][1] * sds[Next(i)][1], -2));
                    denominator += blf.Last();
                }
                for (int i = 0; i < blf.Count; i++)
                {
                    blf[i] = blf[i] / denominator;
                }
            }

            return blf;
        }

        protected int Next(int i, int j = 1) { return (i + j) % n_; }
        protected int Prev(int i, int j = 1) { return (i + n_ - j) % n_; }


    }
}
