using Grasshopper.Kernel.Geometry.Delaunay;
using Rhino.Geometry;
using SCAD;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using System.Windows.Controls.Ribbon;
using System.Windows.Media.Media3D;
using static SCAD.Domain;
using static SCAD.Extensions;
using static SCAD.Utilities;

namespace SCAD
{
    //public class Surface
    //    public class Surface<D, P, R>
    //where D : Domain
    //where P : Parametrization
    //where R : Ribbon     
    //: ISerializable, IEquatable<Point3d>, IComparable<Point3d>, IComparable, IEpsilonComparable<Point3d>, ICloneable, IValidable, IFormattable


    public class Surface<R>
        //where D : Domain, new()
        //where P : Parametrization
        where R : Ribbon, new()
    {

        Domain domain_ = new Domain();
        Parametrization P = new Parametrization();
        //List<Ribbon> R = new List<Ribbon>();

        //protected D domain_;
        //protected P param_;
        //protected List<R> ribbons_;

        protected int n_;
        protected Parametrization param_ = default;
        List<R> ribbons_ = new List<R>();
        private Mesh domainMesh = new Mesh();
        private Brep sphere = new Brep();
        Parametrization_Method pm = new Parametrization_Method();
        Blending_Method bm = new Blending_Method();

        public Surface(Domain_Method dm, Parametrization_Method pm, Blending_Method bm)
        {
            n_ = 0; use_gamma_ = true;
            SetDomain(dm);
            this.pm = pm;
            this.bm = bm;
        }

        private void SetSphere(Brep sphere)
        {
            this.sphere = sphere;
        }

        private void SetDomain(Domain_Method dm)
        {
            if (dm == Domain_Method.Domain_Concave)
                domain_ = new DomainConcave();
            else if (dm == Domain_Method.Domain_Regular)
                domain_ = new DomainRegular();
            else if (dm == Domain_Method.Domain_Curved)
                domain_ = new DomainCurved();
        }

        const double epsilon = 1.0e-8;
        private SurfacePatch sP = new SurfacePatch();
        public Surface() { n_ = 0; use_gamma_ = true; }

        public int N { get { return n_; } }
        ~Surface() { }

        //public Surface(Surface other)       
        //{
        //   //return default;
        //}
        //public static Point3d operator -(Point3d point, Vector3d vector)
        //{
        //    return new Point3d(point.m_x - vector.m_x, point.m_y - vector.m_y, point.m_z - vector.m_z);
        //}

        //public static Surface operator =(Surface srf)
        //{
        //    return new Surface(srf);
        //}

        public Domain GetDomain { get { return domain_; } }

        public Mesh GetDomainMesh { get { return domainMesh; } }

        public Parametrization GetParametrization { get { return P; } }

        public TriMesh MeshTopology(int resolution)
        {
            var mesh = domain_.MeshTopology(resolution);
            return mesh;
        }

        public List<Point2d> Parameters(int resolution)
        {
            var uvs = domain_.Parameters(resolution);
            return uvs;
        }

        public void SetGamma(bool use_gamma)
        {
            use_gamma_ = use_gamma;
        }

        public void SetCurve(int i, NurbsCurve curve)
        {
            if (n_ <= i)
            {
                ribbons_.resize(i + 1);
                n_ = i + 1;
            }
            ribbons_[i] = NewRibbon();
            ribbons_[i].SetCurve(curve);
            domain_.SetSide(i, curve);
        }

        public void SetCurves(List<Curve> curves, List<Curve> tcurves, Brep sphere)
        {
            ribbons_.Clear();
            ribbons_.Capacity = (curves.Count);
            for (int i = 0; i < curves.Count; i++)
            {
                ribbons_.Add(NewRibbon());
                ribbons_.Last().SetCurve(curves[i].ToNurbsCurve());

            }
            //foreach (NurbsCurve curve in curves)
            //{

            //}

            domain_.SetSides(curves);
            domain_.SetRibbonSides(tcurves);
            domain_.SetSphere(sphere);
            SetSphere(sphere);
            n_ = curves.Count;
        }

        public void SetCurves(List<Curve> curves)
        {
            ribbons_.Clear();
            ribbons_.Capacity = (curves.Count);
            for (int i = 0; i < curves.Count; i++)
            {
                ribbons_.Add(NewRibbon());
                ribbons_.Last().SetCurve(curves[i].ToNurbsCurve());

            }
            //foreach (NurbsCurve curve in curves)
            //{

            //}

            domain_.SetSides(curves);
            n_ = curves.Count;
        }

        public virtual void SetupLoop()
        {
            // Tasks:
            // - propagate adjacency information
            // - normalize curves
            // - reverse curves when needed (and normalize once again, for safety)
            for (int i = 0; i < n_; ++i)
                ribbons_[i].Curve().Reparameterize();
            for (int i = 0; i < n_; ++i)
            {
                Ribbon rp = ribbons_[Prev(i)], rn = ribbons_[Next(i)];
                ribbons_[i].SetNeighbors(rp, rn);
                if (i == 0)
                {
                    Point3d r_start = ribbons_[i].Curve().PointAt(0.0);
                    Point3d r_end = ribbons_[i].Curve().PointAt(1.0);
                    Point3d rn_start = rn.Curve().PointAt(0.0);
                    Point3d rn_end = rn.Curve().PointAt(1.0);
                    double end_to_start = (r_end - rn_start).Length;
                    double end_to_end = (r_end - rn_end).Length;
                    double start_to_start = (r_start - rn_start).Length;
                    double start_to_end = (r_start - rn_end).Length;
                    if (Math.Min(start_to_start, start_to_end) < Math.Min(end_to_start, end_to_end))
                    {
                        ribbons_[i].Curve().Reverse();
                        ribbons_[i].Curve().Reparameterize();
                    }
                }
                else
                {
                    Point3d r_start = ribbons_[i].Curve().PointAt(0.0);
                    Point3d r_end = ribbons_[i].Curve().PointAt(1.0);
                    Point3d rp_end = rp.Curve().PointAt(1.0);
                    if ((r_end - rp_end).Length < (r_start - rp_end).Length)
                    {
                        ribbons_[i].Curve().Reverse();
                        ribbons_[i].Curve().Reparameterize();
                    }
                }
            }
        }

        public Vector3d RibbonHandler(int i)
        {
            return ribbons_[i].Handler();
        }

        public void SetRibbonHandler(int i, Vector3d h)
        {
            ribbons_[i].SetHandler(h);
        }

        public void OverrideNormalFence(int i, Vector3d fence)
        {
            ribbons_[i].OverrideNormalFence(fence);
        }

        public double RibbonMultiplier(int i)
        {
            return ribbons_[i].Multiplier();
        }

        public void SetRibbonMultiplier(int i, double m)
        {
            ribbons_[i].SetMultiplier(m);
        }

        public void ResetRibbon(int i)
        {
            ribbons_[i].Reset();
        }
        public virtual void Update(int i)
        {
            if (domain_.Update())
            {
                //param_.update();
            }
            ribbons_[i].Update();
            UpdateCorner(Prev(i));
            UpdateCorner(i);
        }

        public virtual void Update()
        {
            if (domain_.Update())
            {
                //param_.update();
            }
            foreach (var r in ribbons_)
            {
                r.Update();
            }
            UpdateCorners();
        }

        public Domain Domain()
        {
            return domain_;
        }

        public Ribbon Ribbon(int i)
        {
            return ribbons_[i];
        }

        public virtual Point3d Eval(Point2d uv)
        {
            Point3d point_out = sP.GetSurfacePoint(uv.X, uv.Y);
            return point_out;
            //return gen_out;
            //return new Point3d();
        }

        public SurfacePatch GetSurfacePatch { get { return sP; } } 
      




        //public virtual Mesh eval(int resolution)
        //{
        //    TriMesh mesh = domain_.MeshTopology(resolution);
        //    List<Point2d> uvs = domain_.Parameters(resolution);
        //    List<Point3d> points = new List<Point3d>();
        //    points.Capacity = uvs.Count;
        //    foreach (var uv in uvs)
        //    {
        //        points.Add(eval(uv));
        //    }
        //    var msh = mesh.Getmesh;
        //    msh.Vertices.AddVertices(points);
        //    return msh;
        //}

        public virtual Mesh Eval(int resolution)
        {
            TriMesh mesh = domain_.MeshTopology(resolution);
            List<Point2d> uvs = domain_.Parameters(resolution);
            List<Point3d> points = new List<Point3d>();
            points.Capacity = uvs.Count;

            sP = new SurfacePatch(domain_, GetRibbons, corner_data_, pm, bm);

            foreach (var uv in uvs)
            {
                points.Add(Eval(uv));
            }

            domainMesh = mesh.Getmesh;
            var msh = mesh.Getmesh.DuplicateMesh();

            domainMesh.Vertices.AddVertices(uvs.Select(x => new Point3d(x.X, x.Y, 0)));

            msh.Vertices.AddVertices(points);
            P = sP.GetParametrization;
            return msh;
        }

        public virtual Mesh Eval(Mesh rhinoMesh)
        {
            Mesh mesh3d = rhinoMesh.DuplicateMesh();
            var uvs = rhinoMesh.Vertices.ToPoint3dArray().Select(pts3d => new Point2d(pts3d.X, pts3d.Y)).ToList();
            List<Point3d> points = new List<Point3d>();
            points.Capacity = uvs.Count;

            sP = new SurfacePatch(domain_, GetRibbons, corner_data_, pm, bm);
            foreach (var uv in uvs)
            {
                points.Add(Eval(uv));
            }
            for (int i = 0; i < uvs.Count; i++)
            {
                mesh3d.Vertices[i] = new Point3f((float)points[i].X, (float)points[i].Y, (float)points[i].Z);
            }

            P = sP.GetParametrization;
            return mesh3d;
        }
        protected virtual R NewRibbon() { return new R(); }


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
            return ribbons_[i].Eval(new Point2d(si, di));
        }
        protected Point3d CornerInterpolant(int i, List<Point2d> sds)
        {
            double si = sds[i][0], si1 = sds[Next(i)][0];
            var vec = SideInterpolant(i, si, si1) + SideInterpolant(Next(i), si1, 1.0 - si)
              - CornerCorrection(i, 1.0 - si, si1);
            return new Point3d(vec.X, vec.Y, vec.Z);
        }
        protected Point3d CornerInterpolantD(int i, List<Point2d> sds)
        {
            double di = sds[i][1], di1 = sds[Next(i)][1];
            var vec = SideInterpolant(i, 1.0 - di1, di) + SideInterpolant(Next(i), di, di1)
              - CornerCorrection(i, di1, di);
            return new Point3d(vec.X, vec.Y, vec.Z);
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
        protected List<double> BlendSideSingular(List<Point2d> sds)
        {
            List<double> blf = new List<double>(); blf.Capacity = n_;

            int close_to_boundary = 0;
            foreach (var sd in sds)
            {
                if (sd[1] < epsilon)
                    ++close_to_boundary;
            }

            if (close_to_boundary > 0)
            {
                double blend_val = 1.0 / close_to_boundary;
                foreach (var sd in sds)
                {
                    blf.Add(sd[1] < epsilon ? blend_val : 0.0);
                }
            }
            else
            {
                double denominator = 0.0;
                foreach (var sd in sds)
                {
                    blf.Add(Math.Pow(sd[1], -2));
                    denominator += blf.Last();
                }
                for (int i = 0; i < blf.Count; i++)
                {
                    blf[i] = blf[i] / denominator;
                }
            }

            return blf;
        }

        protected List<double> BlendCornerDeficient(List<Point2d> sds)
        {
            List<double> blf = new List<double>(); blf.Capacity = n_;
            for (int i = 0; i < n_; ++i)
            {
                int ip = Next(i);
                if (sds[i][1] < epsilon && sds[ip][1] < epsilon)
                {
                    blf.Add(1.0);
                    continue;
                }
                blf.Add((sds[ip][1] * hermite(0, 1.0 - sds[i][0]) * hermite(0, sds[i][1]) +
                               sds[i][1] * hermite(0, sds[ip][0]) * hermite(0, sds[ip][1])) /
                              (sds[i][1] + sds[ip][1]));
            }
            return blf;
        }

        protected int Next(int i, int j = 1) { return (i + j) % n_; }
        protected int Prev(int i, int j = 1) { return (i + n_ - j) % n_; }




        public List<Ribbon> GetRibbons { get { return ribbons_.Select(x => (Ribbon)x).ToList(); } }


        //private struct CornerData
        //{
        //    public CornerData(Point3d point, Vector3d tangent1, Vector3d tangent2, Vector3d twist1, Vector3d twist2)
        //    {
        //        this.point = point;
        //        this.tangent1 = tangent1;
        //        this.tangent2 = tangent2;
        //        this.twist1 = twist1;
        //        this.twist2 = twist2;
        //    }
        //    public Point3d point;
        //    public Vector3d tangent1, tangent2, twist1, twist2;
        //};

        //private void UpdateCorner(int i)
        //{
        //    var curves = domain_.Curves;
        //    const double step = 1.0e-4;
        //    int ip = Next(i);

        //    sphere.Surfaces[0].ClosestPoint(curves[i].PointAt(1), out double longparam, out double latparam);
        //    Vector3d sphereVec = Vector3d.CrossProduct(sphere.Surfaces[0].NormalAt(longparam, latparam), curves[i].TangentAt(1));
        //    sphereVec.Unitize();

        //    Point3d point;
        //    Vector3d tangent1, tangent2, twist1, twist2;
        //    Vector3d der;
        //    Vector3d d1, d2;
        //    der = sphereVec;
        //    point = curves[i].PointAt(1.0);
        //    tangent1 = der;//-der

        //    sphere.Surfaces[0].ClosestPoint(curves[ip].PointAt(0), out double longparamip, out double latparamip);
        //    Vector3d sphereVecIP = Vector3d.CrossProduct(sphere.Surfaces[0].NormalAt(longparamip, latparamip), curves[ip].TangentAt(0));
        //    sphereVecIP.Unitize();

        //    der = sphereVecIP;
        //    tangent2 = der;

        //    sphere.Surfaces[0].ClosestPoint(curves[i].PointAt(1), out longparam, out latparam);
        //    Vector3d sphereNorm1 = sphere.Surfaces[0].NormalAt(longparam, latparam);
        //    sphere.Surfaces[0].ClosestPoint(curves[i].PointAt(1 - step), out longparam, out latparam);
        //    Vector3d sphereNorm2 = sphere.Surfaces[0].NormalAt(longparam, latparam);

        //    d1 = ribbons_[i].CrossDerivative(1.0, sphereNorm1);
        //    d2 = ribbons_[i].CrossDerivative(1.0 - step, sphereNorm2);
        //    twist1 = (d2 - d1) / step;

        //    sphere.Surfaces[0].ClosestPoint(curves[ip].PointAt(0), out longparam, out latparam);
        //    sphereNorm1 = sphere.Surfaces[0].NormalAt(longparam, latparam);
        //    sphere.Surfaces[0].ClosestPoint(curves[ip].PointAt(step), out longparam, out latparam);
        //    sphereNorm2 = sphere.Surfaces[0].NormalAt(longparam, latparam);

        //    d1 = ribbons_[ip].CrossDerivative(0.0, sphereNorm1);
        //    d2 = ribbons_[ip].CrossDerivative(step, sphereNorm2);
        //    twist2 = (d2 - d1) / step;
        //    corner_data_[i] = new CornerData(point, tangent1, tangent2, twist1, twist2);
        //}


        private void UpdateCorner(int i)
        {
            const double step = 1.0e-4;
            int ip = Next(i);

            Point3d point;
            Vector3d tangent1, tangent2, twist1, twist2;
            Vector3d[] der;
            Vector3d d1, d2;
            CornerData data = new CornerData();
            der = ribbons_[i].Curve().DerivativeAt(1.0, 1);
            point = ribbons_[i].Curve().PointAt(1.0);
            tangent1 = -der[1];
            der = ribbons_[ip].Curve().DerivativeAt(0.0, 1);
            tangent2 = der[1];
            d1 = ribbons_[i].CrossDerivative(1.0);
            d2 = ribbons_[i].CrossDerivative(1.0 - step);
            twist1 = (d2 - d1) / step;
            d1 = ribbons_[ip].CrossDerivative(0.0);
            d2 = ribbons_[ip].CrossDerivative(step);
            twist2 = (d2 - d1) / step;
            corner_data_[i] = new CornerData(point, tangent1, tangent2, twist1, twist2);
        }

        private void UpdateCorners()
        {
            corner_data_.resize(n_);
            for (int i = 0; i < n_; ++i)
                UpdateCorner(i);
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
        private List<CornerData> corner_data_ = new List<CornerData>();
        private bool use_gamma_;
    }
}
