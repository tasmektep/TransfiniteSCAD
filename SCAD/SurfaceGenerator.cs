using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
using Rhino.Geometry.Collections;

namespace SCAD
{

    //public class domain
    //{
    //    int n_;
    //    double M_PI = Math.PI;
    //    const double epsilon = 1.0e-8;
    //    Point2d center_;

    //    List<Curve> curves_;
    //    List<Point2d> parameters_ = new List<Point2d>();
    //    List<Point2d> vertices_ = new List<Point2d>();
    //    List<Vector2d> du_ = new List<Vector2d>(), dv_ = new List<Vector2d>();

    //    public List<Curve> curves { get { return curves_; } }
    //    public domain() { n_ = 0; }

    //    public (int,int) DomainSize { get{ return (n_, n_); } }


    //    public domain(int n) { n_ = n; }

    //    ~domain() { }

    //    public void setSide(int i, Curve curve)
    //    {
    //        if (curves_.Count() <= i)
    //            curves_.Capacity = i + 1;
    //        curves_[i] = curve;
    //    }

    //    public void setSides(List<Curve> curves) { curves_ = curves; }

    //    // Computes everything from vertices
    //    // except for parameters, which are computed only when needed
    //    public bool update()
    //    {
    //        n_ = vertices_.Count();
    //        computeCenter();
    //        parameters_.Clear();
    //        du_.resize(n_); dv_.resize(n_);
    //        for (int i = 0; i < n_; ++i)
    //        {
    //            du_[i] = vertices_[i] - vertices_.Prev(i);
    //            dv_[i] = new Vector2d(du_[i].Y, -du_[i].X);
    //            if ((center_ - vertices_[i]) * dv_[i] < 0)
    //                dv_[i] = -dv_[i];
    //        }
    //        return true;
    //    }

    //    public bool Update()
    //    {
    //        int m = curves_.Count();
    //        if (n_ == m)
    //            return false;

    //        if (m == 4)
    //        {
    //            vertices_ = new List<Point2d>() { new Point2d(1, 1), new Point2d(-1, 1), new Point2d(-1, -1), new Point2d(1, -1) };
    //            return update();
    //        }

    //        double alpha = 2.0 * M_PI / m;
    //        vertices_.resize(m);
    //        for (int i = 0; i < m; ++i)
    //        {
    //            //int j = (i + 1) % curves_.Count;
    //            vertices_[i] = new Point2d(curves_[i % curves_.Count].PointAtStart.X, curves_[i % curves_.Count].PointAtStart.Y);

    //            //vertices_[i] = new Point2d(Math.Cos(alpha * i), Math.Sin(alpha * i));

    //        }

    //        return update();
    //    }

    //    public List<Point2d> vertices() { return vertices_; }

    //    public Point2d toLocal(int i, Vector2d v)
    //    {
    //        return new Point2d(v * du_[i], v * dv_[i]) / du_[i].SquareLength;
    //    }

    //    public Point2d fromLocal(int i, Vector2d v)
    //    {
    //        var res = du_[i] * v.X + dv_[i] * v.Y;
    //        return new Point2d(res.X, res.Y);
    //    }

    //    public bool intersectEdgeWithRay(int i, Point2d p, Vector2d v, Point2d result)
    //    {
    //        Point2d q1 = vertices_.Prev(i), q2 = vertices_[i];
    //        double denom = v.X * (q2[1] - q1[1]) - v.Y * (q2[0] - q1[0]);

    //        if (Math.Abs(denom) < epsilon)
    //            return false;

    //        double t = (v.X * (p[1] - q1[1]) - v.Y * (p[0] - q1[0])) / denom;

    //        if (t < -epsilon || t > 1 + epsilon)
    //            return false;

    //        t = inrange(0.0, t, 1.0);
    //        result = q1 * (1 - t) + q2 * t;

    //        return true;
    //    }

    //    public int size() { return n_; }

    //    public virtual int meshSize(int resolution)
    //    {
    //        if (n_ == 3)
    //            return (resolution + 1) * (resolution + 2) / 2;
    //        if (n_ == 4)
    //            return (resolution + 1) * (resolution + 1);
    //        return 1 + n_ * resolution * (resolution + 1) / 2;
    //    }

    //    public List<Point2d> parameters(int resolution)
    //    {
    //        int size = meshSize(resolution);
    //        if (parameters_.Count() == size)
    //            return parameters_;
    //        parameters_ = parametersImpl(resolution);
    //        return parameters_;
    //    }
    //    protected virtual List<Point2d> parametersImpl(int resolution)
    //    {
    //        List<Point2d> parameters = new List<Point2d>(meshSize(resolution));

    //        if (n_ == 3)
    //        {
    //            for (int j = 0; j <= resolution; ++j)
    //            {
    //                double u = (double)j / resolution;
    //                var p = vertices_[0] * u + vertices_[2] * (1 - u);
    //                var q = vertices_[1] * u + vertices_[2] * (1 - u);
    //                for (int k = 0; k <= j; ++k)
    //                {
    //                    double v = j == 0 ? 1.0 : (double)k / j;
    //                    parameters.Add(p * (1 - v) + q * v);
    //                }
    //            }
    //        }
    //        else if (n_ == 4)
    //        {
    //            for (int j = 0; j <= resolution; ++j)
    //            {
    //                double u = (double)j / resolution;
    //                var p = vertices_[0] * (1 - u) + vertices_[1] * u;
    //                var q = vertices_[3] * (1 - u) + vertices_[2] * u;
    //                for (int k = 0; k <= resolution; ++k)
    //                {
    //                    double v = (double)k / resolution;
    //                    parameters.Add(p * (1 - v) + q * v);
    //                }
    //            }
    //        }
    //        else
    //        { // n_ > 4
    //            parameters.Add(center_);
    //            for (int j = 1; j <= resolution; ++j)
    //            {
    //                double u = (double)j / (double)resolution;
    //                for (int k = 0; k < n_; ++k)
    //                    for (int i = 0; i < j; ++i)
    //                    {
    //                        double v = (double)i / (double)j;
    //                        Point2d ep = vertices_.Prev(k) * (1.0 - v) + vertices_[k] * v;
    //                        Point2d p = center_ * (1.0 - u) + ep * u;
    //                        parameters.Add(p);
    //                    }
    //            }
    //        }
    //        return parameters;
    //    }

    //    public bool onEdge(int resolution, int index)
    //    {
    //        if (n_ == 3)
    //        {
    //            if (index >= meshSize(resolution) - resolution - 1)
    //                return true;
    //            bool issquare(int n_)
    //            {
    //                int root = (int)Math.Round(Math.Sqrt(n_));
    //                return root * root == n_;
    //            };
    //            int n = index * 8 + 1;
    //            return issquare(n) || issquare(n + 8);
    //        }
    //        if (n_ == 4)
    //        {
    //            return index <= resolution || index >= (resolution + 1) * resolution ||
    //              index % (resolution + 1) == 0 || index % (resolution + 1) == resolution;
    //        }
    //        return index >= meshSize(resolution) - n_ * resolution;
    //    }

    //    public virtual TriMesh MeshTopology(int resolution)
    //    {
    //        TriMesh mesh = new TriMesh();
    //        mesh.resizePoints(meshSize(resolution));

    //        if (n_ == 3)
    //        {
    //            int prev = 0, current = 1;
    //            for (int i = 0; i < resolution; ++i)
    //            {
    //                for (int j = 0; j < i; ++j)
    //                {
    //                    mesh.addTriangle(current + j, current + j + 1, prev + j);
    //                    mesh.addTriangle(current + j + 1, prev + j + 1, prev + j);
    //                }
    //                mesh.addTriangle(current + i, current + i + 1, prev + i);
    //                prev = current;
    //                current += i + 2;
    //            }
    //        }
    //        else if (n_ == 4)
    //        {
    //            for (int i = 0; i < resolution; ++i)
    //                for (int j = 0; j < resolution; ++j)
    //                {
    //                    int index = i * (resolution + 1) + j;
    //                    mesh.addTriangle(index, index + resolution + 1, index + 1);
    //                    mesh.addTriangle(index + 1, index + resolution + 1, index + resolution + 2);
    //                }
    //        }
    //        else
    //        { // n_ > 4
    //            int inner_start = 0, outer_vert = 1;
    //            for (int layer = 1; layer <= resolution; ++layer)
    //            {
    //                int inner_vert = inner_start, outer_start = outer_vert;
    //                for (int side = 0; side < n_; ++side)
    //                {
    //                    int vert = 0;
    //                    while (true)
    //                    {
    //                        int next_vert = (side == n_ - 1 && vert == layer - 1) ? outer_start : (outer_vert + 1);
    //                        mesh.addTriangle(inner_vert, outer_vert, next_vert);
    //                        ++outer_vert;
    //                        if (++vert == layer)
    //                            break;
    //                        int inner_next = (side == n_ - 1 && vert == layer - 1) ? inner_start : (inner_vert + 1);
    //                        mesh.addTriangle(inner_vert, next_vert, inner_next);
    //                        inner_vert = inner_next;
    //                    }
    //                }
    //                inner_start = outer_start;
    //            }
    //        }
    //        return mesh;
    //    }

    //    public Point2d center() { return center_; }

    //    public virtual Point2d edgePoint(int i, double s)
    //    {
    //        return vertices_[i] * s + vertices_.Prev(i) * (1.0 - s);
    //    }

    //    public double edgeLength(int i)
    //    {
    //        return (vertices_[i] - vertices_.Prev(i)).Length;
    //    }

    //    public double angle(int i)
    //    {
    //        Vector2d v1 = vertices_.Prev(i) - vertices_[i];
    //        Vector2d v2 = vertices_.Next(i) - vertices_[i];
    //        return Math.Acos(inrange(-1, v1.Length * v2.Length, 1));
    //    }

    //    virtual protected void computeCenter()
    //    {
    //        List<double> lengths = new List<double>(n_);
    //        for (int i = 0; i < n_; ++i)
    //            lengths.Add((vertices_[i] - vertices_.Prev(i)).Length);
    //        center_ = new Point2d(0, 0);
    //        for (int i = 0; i < n_; ++i)
    //            center_ += vertices_[i] * (lengths[i] + lengths.Next(i));
    //        center_ /= lengths.Sum() * 2; ;

    //    }
    //    double inrange(double min, double x, double max)
    //    {
    //        if (x < min)
    //            return min;
    //        if (x > max)
    //            return max;
    //        return x;
    //    }

    //}

}
