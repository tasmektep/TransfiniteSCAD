using Grasshopper.Kernel.Geometry.Delaunay;
using Rhino.Geometry;
using SCAD;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Media.Media3D;
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
    public class Surface<D, P, R>
        where D : Domain, new()
        where P : Parametrization
        where R : Ribbon, new()
    {

        const double epsilon = 1.0e-8;

        public Surface() { n_ = 0; use_gamma_ = true; }

        public int n { get { return n_; } }
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

        public void setGamma(bool use_gamma)
        {
            use_gamma_ = use_gamma;
        }
        public void setCurve(int i, NurbsCurve curve)
        {
            if (n_ <= i)
            {
                ribbons_.resize(i + 1);
                n_ = i + 1;
            }
            ribbons_[i] = newRibbon();
            ribbons_[i].setCurve(curve);
            domain_.SetSide(i, curve);
        }

        public void setCurves(List<Curve> curves)
        {
            ribbons_.Clear();
            ribbons_.Capacity = (curves.Count);
            foreach (NurbsCurve curve in curves)
            {
                ribbons_.Add(newRibbon());
                ribbons_.Last().setCurve(curve);
            }

            domain_.SetSides(curves);
            n_ = curves.Count;
        }

        public virtual void setupLoop()
        {
            // Tasks:
            // - propagate adjacency information
            // - normalize curves
            // - reverse curves when needed (and normalize once again, for safety)
            for (int i = 0; i < n_; ++i)
                ribbons_[i].curve().Reparameterize();
            for (int i = 0; i < n_; ++i)
            {
                Ribbon rp = ribbons_[prev(i)], rn = ribbons_[next(i)];
                ribbons_[i].setNeighbors(rp, rn);
                if (i == 0)
                {
                    Point3d r_start = ribbons_[i].curve().PointAt(0.0);
                    Point3d r_end = ribbons_[i].curve().PointAt(1.0);
                    Point3d rn_start = rn.curve().PointAt(0.0);
                    Point3d rn_end = rn.curve().PointAt(1.0);
                    double end_to_start = (r_end - rn_start).Length;
                    double end_to_end = (r_end - rn_end).Length;
                    double start_to_start = (r_start - rn_start).Length;
                    double start_to_end = (r_start - rn_end).Length;
                    if (Math.Min(start_to_start, start_to_end) < Math.Min(end_to_start, end_to_end))
                    {
                        ribbons_[i].curve().Reverse();
                        ribbons_[i].curve().Reparameterize();
                    }
                }
                else
                {
                    Point3d r_start = ribbons_[i].curve().PointAt(0.0);
                    Point3d r_end = ribbons_[i].curve().PointAt(1.0);
                    Point3d rp_end = rp.curve().PointAt(1.0);
                    if ((r_end - rp_end).Length < (r_start - rp_end).Length)
                    {
                        ribbons_[i].curve().Reverse();
                        ribbons_[i].curve().Reparameterize();
                    }
                }
            }
        }





        public Vector3d ribbonHandler(int i)
        {
            return ribbons_[i].handler();
        }

        public void setRibbonHandler(int i, Vector3d h)
        {
            ribbons_[i].setHandler(h);
        }
        public void overrideNormalFence(int i, Vector3d fence)
        {
            ribbons_[i].overrideNormalFence(fence);
        }

        public double ribbonMultiplier(int i)
        {
            return ribbons_[i].multiplier();
        }

        public void setRibbonMultiplier(int i, double m)
        {
            ribbons_[i].setMultiplier(m);
        }

        public void resetRibbon(int i)
        {
            ribbons_[i].reset();
        }

        public virtual void update(int i)
        {
            if (domain_.update())
            {
                //param_.update();
            }
            ribbons_[i].update();
            updateCorner(prev(i));
            updateCorner(i);
        }
        public virtual void update()
        {
            if (domain_.update())
            {
                //param_.update();
            }
            foreach (var r in ribbons_)
            {
                r.update();
            }
            updateCorners();
        }
        public Domain domain()
        {
            return domain_;
        }
        //public Parameterization parameterization(){
        //  return param_;
        //}

        public Ribbon ribbon(int i)
        {
            return ribbons_[i];
        }

        public virtual Point3d eval(Point2d uv)
        {
            return new Point3d();
        }
        public virtual Mesh eval(int resolution)
        {
            TriMesh mesh = domain_.MeshTopology(resolution);
            List<Point2d> uvs = domain_.Parameters(resolution);
            List<Point3d> points = new List<Point3d>();
            points.Capacity = uvs.Count;
            foreach (var uv in uvs)
            {
                points.Add(eval(uv));
            }
            var msh = mesh.Getmesh;
            msh.Vertices.AddVertices(points);
            return msh;
        }
        protected virtual R newRibbon() { return new R(); }


        protected Point3d cornerCorrection(int i, double s1, double s2)
        {
            // Assumes that both s1 and s2 are 0 at the corner,
            // s1 increases towards corner (i-1), and s2 towards corner (i+1).
            s1 = inrange(0, gamma(s1), 1);
            s2 = inrange(0, gamma(s2), 1);
            return corner_data_[i].point
              + corner_data_[i].tangent1 * s1
              + corner_data_[i].tangent2 * s2
              + rationalTwist(s1, s2, corner_data_[i].twist2, corner_data_[i].twist1) * s1 * s2;
        }
        protected Point3d sideInterpolant(int i, double si, double di)
        {
            si = inrange(0, si, 1);
            di = Math.Max(gamma(di), 0.0);
            return ribbons_[i].eval(new Point2d(si, di));
        }
        protected Point3d cornerInterpolant(int i, List<Point2d> sds)
        {
            double si = sds[i][0], si1 = sds[next(i)][0];
            var vec = sideInterpolant(i, si, si1) + sideInterpolant(next(i), si1, 1.0 - si)
              - cornerCorrection(i, 1.0 - si, si1);
            return new Point3d(vec.X, vec.Y, vec.Z);
        }
        protected Point3d cornerInterpolantD(int i, List<Point2d> sds)
        {
            double di = sds[i][1], di1 = sds[next(i)][1];
            var vec = sideInterpolant(i, 1.0 - di1, di) + sideInterpolant(next(i), di, di1)
              - cornerCorrection(i, di1, di);
            return new Point3d(vec.X, vec.Y, vec.Z);
        }
        protected List<double> blendCorner(List<Point2d> sds)
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
                    int ip = next(i);
                    if (close_to_boundary > 1)
                        blf.Add(sds[i][1] < epsilon && sds[ip][1] < epsilon ? 1.0 : 0.0);
                    else if (sds[i][1] < epsilon)
                    {
                        double tmp = Math.Pow(sds[ip][1], -2);
                        blf.Add(tmp / (tmp + Math.Pow(sds[prev(i)][1], -2)));
                    }
                    else if (sds[ip][1] < epsilon)
                    {
                        double tmp = Math.Pow(sds[i][1], -2);
                        blf.Add(tmp / (tmp + Math.Pow(sds[next(ip)][1], -2)));
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
                    blf.Add(Math.Pow(sds[i][1] * sds[next(i)][1], -2));
                    denominator += blf.Last();
                }
                for (int i = 0; i < blf.Count; i++)
                {
                    blf[i] = blf[i] / denominator;
                }
            }

            return blf;
        }
        protected List<double> blendSideSingular(List<Point2d> sds)
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

        protected List<double> blendCornerDeficient(List<Point2d> sds)
        {
            List<double> blf = new List<double>(); blf.Capacity = n_;
            for (int i = 0; i < n_; ++i)
            {
                int ip = next(i);
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

        protected int next(int i, int j = 1) { return (i + j) % n_; }
        protected int prev(int i, int j = 1) { return (i + n_ - j) % n_; }


        protected int n_;
        protected D domain_ = new D();
        protected P param_ = default;
        protected List<R> ribbons_ = new List<R>();

        public List<Ribbon> GetRibbons { get { return ribbons_.Select(x => (Ribbon)x).ToList(); } }

        //protected D domain_;
        //protected P param_;
        //protected List<R> ribbons_;

        private struct CornerData
        {
            public CornerData(Point3d point, Vector3d tangent1, Vector3d tangent2, Vector3d twist1, Vector3d twist2)
            {
                this.point = point;
                this.tangent1 = tangent1;
                this.tangent2 = tangent2;
                this.twist1 = twist1;
                this.twist2 = twist2;
            }
            public Point3d point;
            public Vector3d tangent1, tangent2, twist1, twist2;
        };

        private void updateCorner(int i)
        {
            const double step = 1.0e-4;
            int ip = next(i);

            Point3d point;
            Vector3d tangent1, tangent2, twist1, twist2;
            Vector3d[] der;
            Vector3d d1, d2;
            CornerData data = new CornerData();
            der = ribbons_[i].curve().DerivativeAt(1.0, 1);
            point = ribbons_[i].curve().PointAt(1.0);
            tangent1 = -der[1];
            der = ribbons_[ip].curve().DerivativeAt(0.0, 1);
            tangent2 = der[1];
            d1 = ribbons_[i].crossDerivative(1.0);
            d2 = ribbons_[i].crossDerivative(1.0 - step);
            twist1 = (d2 - d1) / step;
            d1 = ribbons_[ip].crossDerivative(0.0);
            d2 = ribbons_[ip].crossDerivative(step);
            twist2 = (d2 - d1) / step;
            corner_data_[i] = new CornerData(point, tangent1, tangent2, twist1, twist2);
        }

        private void updateCorners()
        {
            corner_data_.resize(n_);
            for (int i = 0; i < n_; ++i)
                updateCorner(i);
        }
        private double gamma(double d)
        {
            if (use_gamma_)
                return d / (2.0 * d + 1.0);
            return d;
        }
        private static Vector3d rationalTwist(double u, double v, Vector3d f, Vector3d g)
        {
            if (Math.Abs(u + v) < epsilon)
                return new Vector3d(0, 0, 0);
            return (f * u + g * v) / (u + v);
        }
        private List<CornerData> corner_data_ = new List<CornerData>();
        private bool use_gamma_;
    }
}
