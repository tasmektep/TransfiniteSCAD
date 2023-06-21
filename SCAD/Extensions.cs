using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SCAD
{
    public static class Extensions
    {
        public struct CornerData
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


        public static void Reparameterize(this NurbsCurve curve)
        {
            Interval I = new Interval(0, 1);
            curve.Domain = I;

        }
        public static double inrange(double min, double x, double max)
        {
            if (x < min)
                return min;
            if (x > max)
                return max;
            return x;
        }

        public static int IndexWrapper(int index, int list_count)
        {
            return ((index % list_count) + list_count) % list_count;
        }


        public static T Next<T>(this List<T> list, int index)
        {
            return list[(index + 1) % list.Count];
        }

        public static T Prev<T>(this List<T> list, int index)
        {
            return list[(index + list.Count - 1) % list.Count];
        }

        public static NurbsCurve Create(this NurbsCurve @this, IEnumerable<Point3d> points, bool periodic = false, int degree = -1)
        {
            if (degree == -1)
                return NurbsCurve.Create(periodic, points.Count() - 1, points);
            else
                return NurbsCurve.Create(periodic, degree, points);

        }
    }

    public static class Utilities
    {
        public static double hermite(int i, double t)
        {
            switch (i)
            {
                case 0: return Math.Pow(1 - t, 3) + 3.0 * Math.Pow(1 - t, 2) * t;
                case 1: return Math.Pow(1 - t, 2) * t;
                case 2: return (1 - t) * Math.Pow(t, 2);
                case 3: return 3.0 * (1 - t) * Math.Pow(t, 2) + Math.Pow(t, 3);
            }
            return -1.0;                  // should not come here
        }

    }

    public class TriMesh
    {
        List<Point3d> m_vertices;
        Mesh mesh;
        private int n_ = 3;

        public Mesh Getmesh { get { return mesh; } }
        public void resizePoints(int n) { m_vertices = new List<Point3d>(n); }

        public TriMesh() { mesh = new Mesh(); }

        public void setPoints(List<Point3d> pts)
        {
            mesh.Vertices.AddVertices(pts);
        }

        public void addTriangle(int vt_1, int vt_2, int vt_3) { this.mesh.Faces.AddFace(vt_1, vt_2, vt_3); }
        public int meshSize(int resolution)
        {
            if (n_ == 3)
                return (resolution + 1) * (resolution + 2) / 2;
            if (n_ == 4)
                return (resolution + 1) * (resolution + 1);
            return 1 + n_ * resolution * (resolution + 1) / 2;
        }

    }

    //public static class PointÊxtension
    //{
    //    public static List<Point3d> to3d(this List<Point3d> ds)
    //    {
    //        List<Point3d> vertices3d = new List<Point3d>();
    //        foreach (var item in ds)
    //        {
    //            vertices3d.Add(new Point3d(item.X, item.Y, 0));
    //        }
    //        return vertices3d;
    //    }

    //}



    public static class ListExtension
    {
        private static void resize<T>(this List<T> list, int sz, T c)
        {
            int cur = list.Count;
            if (sz < cur)
                list.RemoveRange(sz, cur - sz);
            else if (sz > cur)
            {
                if (sz > list.Capacity)//this bit is purely an optimisation, to avoid multiple automatic capacity changes.
                    list.Capacity = sz;
                list.AddRange(Enumerable.Repeat(c, sz - cur));
            }
        }
        public static void resize<T>(this List<T> list, int sz) where T : new()
        {
            resize(list, sz, new T());
        }
    }
}


