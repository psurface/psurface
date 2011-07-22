/**
 * @file Box.hh
 * @brief
 */

#ifndef BOX_H
#define BOX_H

#include <algorithm>

// Check for VC9 / VS2008 with installed feature pack.
#if defined(_MSC_VER) && (_MSC_VER>=1500)
    #if defined(_CPPLIB_VER) && _CPPLIB_VER>=505
        #include <array>
    #else
        #error Please install the Visual C++ 2008 Feature Pack for TR1 support.
    #endif
#else
    #include <tr1/array>
#endif

#ifndef PSURFACE_STANDALONE
#include <mclib/McVec3f.h>
#endif

/** \brief A axis-parallel box in a Euclidean space
    \tparam C Type used for coordinate components
    \tparam dim Dimension of the box
*/
template<typename C, int dim>
class Box
{
public:

    /** \brief Default constructor.  Box is not initialized! */
    Box()
    {}

    /** \brief Set box from two points */
    Box(const std::tr1::array<C,dim>& lower, const std::tr1::array<C,dim>& upper) 
        : _lower(lower), _upper(upper)
    {
        for (int i=0; i<dim; i++) {
            _lower[i] = std::min(lower[i],upper[i]);
            _upper[i] = std::max(lower[i],upper[i]);
        }
    }

#ifndef PSURFACE_STANDALONE
    /** \brief Set box from two McVec3f */
    Box(const McVec3f& lower, const McVec3f& upper) 
    {
        for (int i=0; i<dim; i++) {
            _lower[i] = std::min(lower[i],upper[i]);
            _upper[i] = std::max(lower[i],upper[i]);
        }
    }
#endif

    /** \brief Copy constructor */
    Box(const Box& b) : _lower(b._lower), _upper(b._upper)
    {}

    /** \brief Set up box from to diagonal corners */
    void set(const std::tr1::array<C,dim>& lower, const std::tr1::array<C,dim>& upper)
    {
        for (int i=0; i<dim; i++) {
            _lower[i] = std::min(lower[i],upper[i]);
            _upper[i] = std::max(lower[i],upper[i]);
        }
    }

    /** \brief Test whether box contains a given point */
	bool contains(const std::tr1::array<C,dim>& c) const
	{
		for (int i = 0; i < dim; ++i)
			if (c[i] < this->_lower[i] || c[i] >= this->_upper[i])
				return false;
		return true;
	}

#ifndef PSURFACE_STANDALONE
    /** \brief Test whether box contains a given point */
    bool contains(const McVec3f& c) const
    {
        for (int i = 0; i < dim; ++i)
            if (c[i] < this->_lower[i] || c[i] >= this->_upper[i])
                return false;
        return true;
    }
#endif

	bool intersects(const Box& b)
	{
		for (int i = 0; i < dim; ++i)
			if (this->_lower[i] >= b._upper[i] || b._lower[i] >= this->_upper[i])
				return false;
		return true;
	}

    /// Returns intersection of two boxes.
    Box<C,dim> intersectWith(const Box<C,dim> &other) const {

        std::tr1::array<C,dim> zero;
        zero.assign(0);

        Box<C,dim> innerBox(zero,zero);

        for (int i = 0; i < dim; i++) {

            if ((upper()[i] < other.lower()[i]) || (lower()[i] > other.upper()[i]))
                return Box<C,dim>(zero,zero);
            
            innerBox._lower[i] = std::max(lower()[i],other.lower()[i]);
            innerBox._upper[i] = std::min(upper()[i],other.upper()[i]);
        }

        return innerBox;
    }

	std::tr1::array<C,dim> center() const
	{
            std::tr1::array<C,dim> center;
            for (int i = 0; i < dim; ++i)
                center[i] = 0.5*(_upper[i]+_lower[i]);
            return center;
	}


	double size(int i)
	{
		return _upper[i]-_lower[i];
	}

        std::tr1::array<C,dim>& lower()
	{
		return _lower;
	}

	std::tr1::array<C,dim>& upper()
	{
		return _upper;
	}

	const std::tr1::array<C,dim>& lower() const
	{
		return _lower;
	}

	const std::tr1::array<C,dim>& upper() const
	{
		return _upper;
	}

        /// Extends the box to contain given point.
    void extendBy(const std::tr1::array<C,dim>& point){
        for (int i=0; i<dim; i++) {
            _lower[i] = std::min(_lower[i], point[i]);
            _upper[i] = std::max(_upper[i], point[i]);
        }
    }

#ifndef PSURFACE_STANDALONE
    void extendBy(const McVec3f& point){
        for (int i=0; i<dim; i++) {
            _lower[i] = std::min(_lower[i], point[i]);
            _upper[i] = std::max(_upper[i], point[i]);
        }
    }
#endif

    /// Enlarges the box by 'eps' to each side
    void extendByEps(float eps){
        for (int i=0; i<dim; i++) {
            _lower[i] -= eps;
            _upper[i] += eps;
        }
    }


private:

	std::tr1::array<C,dim> _lower;
	std::tr1::array<C,dim> _upper;
};

#endif // BOX_HH_
