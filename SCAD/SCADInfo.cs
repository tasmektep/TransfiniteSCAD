using Grasshopper;
using Grasshopper.Kernel;
using System;
using System.Drawing;

namespace SCAD
{
    public class SCADInfo : GH_AssemblyInfo
    {
        public override string Name => "Transfinite";

        //Return a 24x24 pixel bitmap to represent this GHA library.
        public override Bitmap Icon => null;

        //Return a short string describing the purpose of this GHA library.
        public override string Description => "";

        public override Guid Id => new Guid("448CAB9E-CD5D-4119-AC5F-6E3DC55A6448");

        //Return a string identifying you or your company.
        public override string AuthorName => "";

        //Return a string representing your preferred contact details.
        public override string AuthorContact => "";
    }
}