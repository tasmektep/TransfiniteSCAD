using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;

namespace SCAD
{
    public class Domain
    {

        int n_;
        double M_PI = Math.PI;
        const double epsilon = 1.0e-8;
        Point2d center_;

        List<Curve> curves_;
        int resolution;
        List<Point2d> parameters_ = new List<Point2d>();
        List<Point2d> vertices_ = new List<Point2d>();
        List<Vector2d> du_ = new List<Vector2d>(), dv_ = new List<Vector2d>();

        private List<Curve> m_BndCurves;
        private List<Curve> m_DomainCurves;

        /// <summary>
        /// Gets boundary curves to create Domain Curves
        /// returns: The Domain Curves
        /// </summary>
        /// <param name="BndCurves"> Boundary Curves </param>
        /// <param name="DomainCurves"> Domain Curves </param>
        public Domain(List<Curve> BndCurves, out List<Curve> DomainCurves)
        {
            m_BndCurves = BndCurves;
            ComputeDomainPolygon();
            DomainCurves = m_DomainCurves;
        }

        private void ComputeDomainPolygon()
        {
            /*double TotL = 0, L;
            List<double> Ls = new List<double>();
            for (int i = 0; i < m_BndCurves.Count; i++)
            {
                L = m_BndCurves[i].GetLength();
                TotL = TotL + L;
                Ls.Add(L);
            }


            List<double> Alphas = new List<double>();
            Alphas.Add(0.0);
            for (int i = 1; i < m_BndCurves.Count; i++)
            {
                L = 0;
                for (int j = 0; j <= i - 1; j++)
                    L = L + Ls[j];
                Alphas.Add(2 * Math.PI * L / TotL);
            }

            m_DomainCurves = new List<Curve>();
            Point3d Pt1, Pt2;
            int k;
            Line Ln;
            for (int i = 0; i < Alphas.Count; i++)
            {
                k = (i + 1) % Alphas.Count;
                Pt1 = new Point3d(Math.Cos(Alphas[i]), Math.Sin(Alphas[i]), 0);
                Pt2 = new Point3d(Math.Cos(Alphas[k]), Math.Sin(Alphas[k]), 0);

                Ln = new Line(Pt1, Pt2);
                m_DomainCurves.Add(Ln.ToNurbsCurve());
            }*/

            Point3d Pt1, Pt2;
            int j;
            Line Ln;
            for (int i = 0; i < m_BndCurves.Count; i++)
            {
                j = (i + 1) % m_BndCurves.Count;

                Pt1 = new Point3d(m_BndCurves[i].PointAtStart.X, m_BndCurves[i].PointAtStart.Y, 0);
                Pt2 = new Point3d(m_BndCurves[j].PointAtStart.X, m_BndCurves[j].PointAtStart.Y, 0);

                Ln = new Line(Pt1, Pt2);
                m_DomainCurves.Add(Ln.ToNurbsCurve());
            }
        }

        //private

        public (Interval bx, Interval by) Bounds
        {
            get
            {
                var multiplayer = 1;
                //Interval x = new Interval(vertices_.Min(pt => pt.X), vertices_.Max(pt => pt.X));
                //Interval y = new Interval(vertices_.Min(pt => pt.Y), vertices_.Max(pt => pt.Y));
                Interval x = new Interval(Math.Floor(vertices_.Min(pt => pt.X)) * multiplayer, Math.Ceiling(vertices_.Max(pt => pt.X)) * multiplayer);
                Interval y = new Interval(Math.Floor(vertices_.Min(pt => pt.Y)) * multiplayer, Math.Ceiling(vertices_.Max(pt => pt.Y)) * multiplayer);
                return (x, y);
            }
        }

        public List<Curve> Curves { get { return curves_; } }

        //public (int, int) DomainSize { get { return (n_, n_); } }

        public List<Point2d> Vertices { get { return vertices_; } }
        public List<Point3d> Vertices3d
        {
            get
            {
                List<Point3d> vertices3d = new List<Point3d>();
                foreach (var item in vertices_)
                {
                    vertices3d.Add(new Point3d(item.X, item.Y, 0));
                }
                return vertices3d;
            }
        }

        public int Resulation { get { return resolution; } }


        public Point2d Center { get { return center_; } }

        public int Size { get { return n_; } }


        //public Domain(int n) { n_ = n; }

        public Domain(int resulation) { this.resolution = resulation; }

        ~Domain() { }

        public void SetSide(int i, Curve curve)
        {
            if (curves_.Count() <= i)
                curves_.Capacity = i + 1;
            curves_[i] = curve;
        }

        public void SetSides(List<Curve> curves) { curves_ = curves; }

        // Computes everything from vertices
        // except for parameters, which are computed only when needed
        private bool update()
        {
            n_ = vertices_.Count();
            ComputeCenter();
            parameters_.Clear();
            du_.resize(n_); dv_.resize(n_);
            for (int i = 0; i < n_; ++i)
            {
                du_[i] = vertices_[i] - vertices_.Prev(i);
                dv_[i] = new Vector2d(du_[i].Y, -du_[i].X);
                if ((center_ - vertices_[i]) * dv_[i] < 0)
                    dv_[i] = -dv_[i];
            }
            return true;
        }

        public bool Update()
        {
            int m = curves_.Count();
            if (n_ == m)
                return false;

            if (m == 4)
            {
                vertices_ = new List<Point2d>() { new Point2d(1, 1), new Point2d(-1, 1), new Point2d(-1, -1), new Point2d(1, -1) };
                return update();
            }

            double alpha = 2.0 * M_PI / m;
            vertices_.resize(m);
            for (int i = 0; i < m; ++i)
            {
                //int j = (i + 1) % curves_.Count;
                //vertices_[i] = new Point2d(Math.Cos(alpha * i), Math.Sin(alpha * i));
                
                vertices_[i] = new Point2d(curves_[i % curves_.Count].PointAtStart.X, curves_[i % curves_.Count].PointAtStart.Y);

            }

            return update();
        }


        public Point2d ToLocal(int i, Vector2d v)
        {
            return new Point2d(v * du_[i], v * dv_[i]) / du_[i].SquareLength;
        }

        public Point2d FromLocal(int i, Vector2d v)
        {
            var res = du_[i] * v.X + dv_[i] * v.Y;
            return new Point2d(res.X, res.Y);
        }

        public bool IntersectEdgeWithRay(int i, Point2d p, Vector2d v, Point2d result)
        {
            Point2d q1 = vertices_.Prev(i), q2 = vertices_[i];
            double denom = v.X * (q2[1] - q1[1]) - v.Y * (q2[0] - q1[0]);

            if (Math.Abs(denom) < epsilon)
                return false;

            double t = (v.X * (p[1] - q1[1]) - v.Y * (p[0] - q1[0])) / denom;

            if (t < -epsilon || t > 1 + epsilon)
                return false;

            t = Inrange(0.0, t, 1.0);
            result = q1 * (1 - t) + q2 * t;

            return true;
        }


        public virtual int MeshSize()
        {
            if (n_ == 3)
                return (resolution + 1) * (resolution + 2) / 2;
            if (n_ == 4)
                return (resolution + 1) * (resolution + 1);
            return 1 + n_ * resolution * (resolution + 1) / 2;
        }

        public List<Point2d> Parameters()
        {
            int size = MeshSize();
            if (parameters_.Count() == size)
                return parameters_;
            parameters_ = ParametersImpl();
            return parameters_;
        }
        protected virtual List<Point2d> ParametersImpl()
        {
            List<Point2d> parameters = new List<Point2d>(MeshSize());

            if (n_ == 3)
            {
                for (int j = 0; j <= resolution; ++j)
                {
                    double u = (double)j / resolution;
                    var p = vertices_[0] * u + vertices_[2] * (1 - u);
                    var q = vertices_[1] * u + vertices_[2] * (1 - u);
                    for (int k = 0; k <= j; ++k)
                    {
                        double v = j == 0 ? 1.0 : (double)k / j;
                        parameters.Add(p * (1 - v) + q * v);
                    }
                }
            }
            else if (n_ == 4)
            {
                for (int j = 0; j <= resolution; ++j)
                {
                    double u = (double)j / resolution;
                    var p = vertices_[0] * (1 - u) + vertices_[1] * u;
                    var q = vertices_[3] * (1 - u) + vertices_[2] * u;
                    for (int k = 0; k <= resolution; ++k)
                    {
                        double v = (double)k / resolution;
                        parameters.Add(p * (1 - v) + q * v);
                    }
                }
            }
            else
            { // n_ > 4
                parameters.Add(center_);
                for (int j = 1; j <= resolution; ++j)
                {
                    double u = (double)j / (double)resolution;
                    for (int k = 0; k < n_; ++k)
                        for (int i = 0; i < j; ++i)
                        {
                            double v = (double)i / (double)j;
                            Point2d ep = vertices_.Prev(k) * (1.0 - v) + vertices_[k] * v;
                            Point2d p = center_ * (1.0 - u) + ep * u;
                            parameters.Add(p);
                        }
                }
            }
            return parameters;
        }

        public bool OnEdge(int index)
        {
            if (n_ == 3)
            {
                if (index >= MeshSize() - resolution - 1)
                    return true;
                bool issquare(int n_)
                {
                    int root = (int)Math.Round(Math.Sqrt(n_));
                    return root * root == n_;
                };
                int n = index * 8 + 1;
                return issquare(n) || issquare(n + 8);
            }
            if (n_ == 4)
            {
                return index <= resolution || index >= (resolution + 1) * resolution ||
                  index % (resolution + 1) == 0 || index % (resolution + 1) == resolution;
            }
            return index >= MeshSize() - n_ * resolution;
        }

        public virtual TriMesh MeshTopology()
        {
            TriMesh mesh = new TriMesh();
            mesh.resizePoints(MeshSize());

            if (n_ == 3)
            {
                int prev = 0, current = 1;
                for (int i = 0; i < resolution; ++i)
                {
                    for (int j = 0; j < i; ++j)
                    {
                        mesh.addTriangle(current + j, current + j + 1, prev + j);
                        mesh.addTriangle(current + j + 1, prev + j + 1, prev + j);
                    }
                    mesh.addTriangle(current + i, current + i + 1, prev + i);
                    prev = current;
                    current += i + 2;
                }
            }
            else if (n_ == 4)
            {
                for (int i = 0; i < resolution; ++i)
                    for (int j = 0; j < resolution; ++j)
                    {
                        int index = i * (resolution + 1) + j;
                        mesh.addTriangle(index, index + resolution + 1, index + 1);
                        mesh.addTriangle(index + 1, index + resolution + 1, index + resolution + 2);
                    }
            }
            else
            { // n_ > 4
                int inner_start = 0, outer_vert = 1;
                for (int layer = 1; layer <= resolution; ++layer)
                {
                    int inner_vert = inner_start, outer_start = outer_vert;
                    for (int side = 0; side < n_; ++side)
                    {
                        int vert = 0;
                        while (true)
                        {
                            int next_vert = (side == n_ - 1 && vert == layer - 1) ? outer_start : (outer_vert + 1);
                            mesh.addTriangle(inner_vert, outer_vert, next_vert);
                            ++outer_vert;
                            if (++vert == layer)
                                break;
                            int inner_next = (side == n_ - 1 && vert == layer - 1) ? inner_start : (inner_vert + 1);
                            mesh.addTriangle(inner_vert, next_vert, inner_next);
                            inner_vert = inner_next;
                        }
                    }
                    inner_start = outer_start;
                }
            }
            return mesh;
        }


        public virtual Point2d edgePoint(int i, double s)
        {
            return vertices_[i] * s + vertices_.Prev(i) * (1.0 - s);
        }

        public double edgeLength(int i)
        {
            return (vertices_[i] - vertices_.Prev(i)).Length;
        }

        public double angle(int i)
        {
            Vector2d v1 = vertices_.Prev(i) - vertices_[i];
            Vector2d v2 = vertices_.Next(i) - vertices_[i];
            return Math.Acos(Inrange(-1, v1.Length * v2.Length, 1));
        }

        virtual protected void ComputeCenter()
        {
            List<double> lengths = new List<double>(n_);
            for (int i = 0; i < n_; ++i)
                lengths.Add((vertices_[i] - vertices_.Prev(i)).Length);
            center_ = new Point2d(0, 0);
            for (int i = 0; i < n_; ++i)
                center_ += vertices_[i] * (lengths[i] + lengths.Next(i));
            center_ /= lengths.Sum() * 2; ;

        }
        double Inrange(double min, double x, double max)
        {
            if (x < min)
                return min;
            if (x > max)
                return max;
            return x;
        }

    }

}
