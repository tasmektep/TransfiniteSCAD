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
    }

    public class TriMesh
    {
        List<Point3d> m_vertices;
        Mesh mesh;
        private int n_ = 3;

        public Mesh Getmesh { get { return mesh; } }
        public void resizePoints(int n) { m_vertices = new List<Point3d>(n); }

        public TriMesh() { mesh = new Mesh(); }

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
        public static void resize<T>(this List<T> list, int sz, T c)
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
