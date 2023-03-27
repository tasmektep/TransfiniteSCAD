using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino.Geometry;
namespace SCAD
{
    class RMF
    {
        const double epsilon = 1.0e-8;

        public struct Frame
        {
            public Frame(double u, double s, Point3d p, Vector3d d, Vector3d n)
            {
                this.d = d;
                this.u = u;
                this.s = s;
                this.p = p;
                this.n = n;
            }
            public double u, s;
            public Point3d p;
            public Vector3d d, n;
        };

        List<Frame> frames_ = new List<Frame>();
        double angleCorrection_;
        Vector3d start_ = new Vector3d();
        Vector3d end_ = new Vector3d();
        const int resolution_ = 100;

        NurbsCurve curve_ = new NurbsCurve(3, 4);
        public void setCurve(NurbsCurve c)
        {
            curve_ = c;
        }

        public void setStart(Vector3d start)
        {
            start_ = start;
        }

        public void setEnd(Vector3d end)
        {
            end_ = end;
        }

        public void update()
        {
            // As described in `Computation of Rotation Minimizing Frames', Wang et al., 2008.
            // A limitation of this method is that it is determined by the starting frame
            // and the curve tangents, so an end frame cannot be supplied.
            frames_.Clear();
            //frames_.reserve(resolution_ + 1);
            List<Vector3d> der; der = curve_.DerivativeAt(0.0, 1).ToList();
            der[1].Unitize();
            Frame f = new Frame(0.0, 0.0, (Point3d)der[0], der[1], start_);
            for (int i = 1; i <= resolution_; ++i)
            {
                double u = (double)i / (double)resolution_;
                frames_.Add(f);
                f = nextFrame(f, u);
            }
            frames_.Add(f);

            // As a workaround, we can add a rotation gradually,
            // by minimizing the total squared angular speed, see section 6.3 in the paper.
            Vector3d rmfEnd = f.n;
            angleCorrection_ = Math.Acos(inrange(-1, end_ * rmfEnd, 1));
            if ((Vector3d.CrossProduct((rmfEnd - end_), end_) * f.d < 0.0))
                angleCorrection_ *= -1.0;
            angleCorrection_ /= f.s;
        }

        public Vector3d eval(double u)
        {
            int i = frames_.FindIndex(x => x.u >= u);
            if (u < frames_[i].u)
                return frames_[i].n;
            Frame f = nextFrame(frames_[i - 1], u);
            rotateFrame(f, f.s * angleCorrection_);
            return f.n;
        }
        double inrange(double min, double x, double max)
        {
            if (x < min)
                return min;
            if (x > max)
                return max;
            return x;
        }
        Matrix rotationMatrix(Vector3d u, double theta)
        {
            Matrix m = new Matrix(3, 3);
            double x = u[0], y = u[1], z = u[2];
            double c = Math.Cos(theta), c1 = 1.0 - c, s = Math.Sin(theta);
            m[0, 0] = c + x * x * c1;
            m[1, 0] = y * x * c1 + z * s;
            m[2, 0] = z * x * c1 - y * s;
            m[0, 1] = x * y * c1 - z * s;
            m[1, 1] = c + y * y * c1;
            m[2, 1] = z * y * c1 + x * s;
            m[0, 2] = x * z * c1 + y * s;
            m[1, 2] = y * z * c1 - x * s;
            m[2, 2] = c + z * z * c1;

            return m;
        }

        void
        rotateFrame(Frame f, double angle)
        {
            Matrix r = rotationMatrix(f.d, angle);

            Vector3d n = new Vector3d(0.0, 0.0, 0.0);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    n[i] += r[i, j] * f.n[j];
            f.n = n;
        }

        Frame nextFrame(Frame prev, double u)
        {
            List<Vector3d> der; der = curve_.DerivativeAt(u, 1).ToList();
            //der[1].normalize();
            Vector3d v1 = der[0] - ((Vector3d)prev.p);
            double c1 = v1 * v1;
            if (c1 < epsilon)
                return prev;
            Vector3d v2 = v1 * 2 / c1;
            Vector3d nL = prev.n - v2 * (v1 * prev.n);
            Vector3d dL = prev.d - v2 * (v1 * prev.d);
            v2 = der[1] - dL;
            double c2 = v2 * v2;
            Vector3d nNext;
            if (c2 < epsilon)
                nNext = nL;
            else
                nNext = nL - v2 * 2 / c2 * (v2 * nL);

            return new Frame(u, prev.s + curve_.GetLength(new Interval(prev.u, u)), (Point3d)der[0], der[1], nNext); ;
        }

    } // namespace Transfinite

}

