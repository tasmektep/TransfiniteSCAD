using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace SCAD
{
    public class Test : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Test class.
        /// </summary>
        public Test()
          : base("Test", "Test",
              "Description",
             "SCAD", "Transfinite")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("", "", "", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("", "", "", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Curve> curves_ = new List<Curve>();
            DA.GetDataList(0, curves_);
            int resolution = 30;
            var dm = new domain();
            dm.setSides(curves_);
            dm.Update();



            var mesh = dm.MeshTopology(resolution);
            var uvs = dm.parameters(resolution);
            DA.SetData(0, mesh.Getmesh);
         
            var msh = mesh.Getmesh;
            var pts = new List<Point3d>();
            var sP = new SurfacePatch();
            for (int i = 0; i < uvs.Count; i++)
            {
                //pts.Add(new Point3d(uvs[i].X, uvs[i].Y, 0));
                Point3d katoout = sP.Kato_Suv(uvs[i].X, uvs[i].Y, curves_);
                pts.Add(katoout);
            }

            msh.Vertices.AddVertices(pts);

        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("5587F090-35B5-45DE-8AE8-64D81549853E"); }
        }
    }
}