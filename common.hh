// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_COMMON_HH
#define DUNE_GRID_IO_FILE_VTK_COMMON_HH

#include <limits>
#include <sstream>
#include <string>


/** @file
    @author Peter Bastian, Christian Engwer
    @brief Common stuff for the VTKWriter

    This file contains common stuff for all instances of VTKWriter.
 */

namespace psurface 
{
  //! \addtogroup VTK
  //! \{

  namespace VTK {

    //////////////////////////////////////////////////////////////////////
    //
    //  VTKOptions
    //

    //! How the bulk data should be stored in the file
    /**
     * \code
     * #include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     */
    enum OutputType {
      //! Output to the file is in ascii.
      ascii,
      //! Output to the file is inline base64 binary.
      base64,
      //! Ouput is to the file is appended raw binary
      appendedraw,
      //! Ouput is to the file is appended base64 binary
      appendedbase64
      // //! Output to the file is compressed inline binary.
      // binarycompressed,
      // //! Ouput is compressed and appended to the file.
      // compressedappended
    };
    //! Whether to produce conforming or non-conforming output.
    /**
     * \code
     * #include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     *
     * This applies to the conformity of the data; a non-conforming grid can
     * still be written in conforming data mode, and it is quite possible for
     * data to be non-conforming on a conforming grid.
     */
    enum DataMode {
      //! Output conforming data.
      /**
       * Neighboring elements share common vertices and thus have a common DoF
       * on that vertex.
       */
      conforming,
      //! Output non-conforming data.
      /**
       * Each element has it's own set of vertices.  The position of a vertex
       * of one element will concide with the position of the corresponding
       * vertex on another element.  This allows for multiple DoFs (one per
       * element) on the "same" vertex.
       */
      nonconforming
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  PrintType
    //

    //! determine a type to safely put another type into a stream
    /**
     * This is mainly interating for character types which should print as
     * their integral value, not as a character.
     */
    template<typename T>
    struct PrintType {
      //! type to convert T to before putting it into a stream with <<
      typedef T Type;
    };

    template<>
    struct PrintType<unsigned char> {
      typedef unsigned Type;
    };

    template<>
    struct PrintType<signed char> {
      typedef int Type;
    };

    template<>
    struct PrintType<char> {
      typedef int Type;
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  TypeName
    //

    //! map type to its VTK name in data array
    /**
     * \tparam T The type whose VTK name is requested
     */
    template<typename T>
    class TypeName {
      static std::string getString() {
        static const unsigned int_sizes[] = { 8, 16, 32, 64, 0 };
        static const unsigned float_sizes[] = { 32, 64, 0 };
        const unsigned* sizes;

        std::ostringstream s;
        if(std::numeric_limits<T>::is_integer) {
          if(std::numeric_limits<T>::is_signed)
            s << "Int";
          else
            s << "UInt";
          sizes = int_sizes;
        }
        else {
          // assume float
          s << "Float";
          sizes = float_sizes;
        }

        static const unsigned size = 8*sizeof(T);
        while(*sizes != 0 && *sizes <= size) ++sizes;
        --sizes;
        s << *sizes;

        return s.str();
      }

    public:
      //! return VTK name of the type
      /**
       * If the type is not known to VTK, return empty string.
       */
      const std::string& operator()() const {
        static const std::string s = getString();
        return s;
      }
    };

    //////////////////////////////////////////////////////////////////////
    //
    //  VTK::GeometryType related stuff
    //

    //! Type representing VTK's entity geometry types
    /**
     * \code
     * #include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     *
     * Only the types which have a corresponding Dune::GeometryType have been
     * included here.  Dune-type names have been used, this mainly makes a
     * difference for vtkPrism, which is known by VTK as VTK_WEDGE.
     */
/*    enum GeometryType {
      vertex = 1,
      line = 3,
      triangle = 5,
      quadrilateral = 9,
      tetrahedron = 10,
      hexahedron = 12,
      prism = 13,
      pyramid = 14
    };
*/
    template<typename T>
    int renumber(const T& t, int i)
    {
      return renumber(t.type(), i);
    }

    //////////////////////////////////////////////////////////////////////
    //
    //  Determine Endianness
    //

    //! determine endianness of this C++ implementation
    /**
     * \returns A string suitable for the byte_order property in VTK files;
     *          either "BigEndian" or "LittleEndian".
     */
    inline std::string getEndiannessString()
    {
      short i = 1;
      if (reinterpret_cast<char*>(&i)[1] == 1)
        return "BigEndian";
      else
        return "LittleEndian";
    }

    //////////////////////////////////////////////////////////////////////
    //
    //  which type of vtkfile to write
    //

    //! which type of VTK file to write
    /**
     * \code
     * #include <dune/grid/io/file/vtk/common.hh>
     * \endcode
     */
    enum FileType {
      //! for .vtp files (PolyData)
      polyData,
      //! for .vtu files (UnstructuredGrid)
      unstructuredGrid
    };

  } // namespace VTK

  //! \} group VTK

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_COMMON_HH
