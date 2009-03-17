/*
 *  Filename:    Box.hh
 *  Version:     1.0
 *  Created on:  Jan 13, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     liboctree
 *  Description: <short_description>
 *  subversion:  $Id$
 *
 */
/**
 * @file Box.hh
 * @brief
 */

#ifndef BOX_H
#define BOX_H

#include <algorithm>
#include <mclib/McVec3f.h>

/** \brief A axis-parallel box in a Euclidean space
    \tparam C Type used for coordinates
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
    Box(const C& lower, const C& upper) : _lower(lower), _upper(upper)
    {
        for (int i=0; i<3; i++) {
            _lower[i] = std::min(lower[i],upper[i]);
            _upper[i] = std::max(lower[i],upper[i]);
        }
    }

    /** \brief Set box from two points */
    Box(const McVec3f& lower, const McVec3f& upper)
    {
        for (int i=0; i<3; i++) {
            _lower[i] = std::min(lower[i],upper[i]);
            _upper[i] = std::max(lower[i],upper[i]);
        }
    }

    /** \brief Copy constructor */
	Box(const Box& b) : _lower(b._lower), _upper(b._upper)
    {}

    /** \brief Set up box from to diagonal corners */
    void set(const C& lower, const C& upper)
    {
        for (int i=0; i<3; i++) {
            _lower[i] = std::min(lower[i],upper[i]);
            _upper[i] = std::max(lower[i],upper[i]);
        }
    }

    /** \brief Set up box from to diagonal corners */
    void set(const McVec3f& lower, const McVec3f& upper)
    {
        for (int i=0; i<3; i++) {
            _lower[i] = std::min(lower[i],upper[i]);
            _upper[i] = std::max(lower[i],upper[i]);
        }
    }

    /** \brief Test whether box contains a given point */
	bool contains(const C& c) const
	{
		for (int i = 0; i < dim; ++i)
			if (c[i] < this->_lower[i] || c[i] >= this->_upper[i])
				return false;
		return true;
	}

    /** \brief Test whether box contains a given point */
    bool contains(const McVec3f& c) const
    {
        for (int i = 0; i < dim; ++i)
            if (c[i] < this->_lower[i] || c[i] >= this->_upper[i])
                return false;
        return true;
    }

	bool intersects(const Box& b)
	{
		for (int i = 0; i < dim; ++i)
			if (this->_lower[i] >= b._upper[i] || b._lower[i] >= this->_upper[i])
				return false;
		return true;
	}

    /// Returns intersection of two boxes.
    Box<C,dim> intersectWith(const Box<C,dim> &other) const {

        C zero;
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

	C center() const
	{
            C center;
            for (int i = 0; i < dim; ++i)
                center[i] = 0.5*(_upper[i]+_lower[i]);
            return center;
	}


	double size(int i)
	{
		return _upper[i]-_lower[i];
	}

	const C& lower() const
	{
		return _lower;
	}

	const C& upper() const
	{
		return _upper;
	}

        /// Extends the box to contain given point.
    void extendBy(const C& point){
        for (int i=0; i<dim; i++) {
            _lower[i] = std::min(_lower[i], point[i]);
            _upper[i] = std::max(_upper[i], point[i]);
        }
    }

    void extendBy(const McVec3f& point){
        for (int i=0; i<dim; i++) {
            _lower[i] = std::min(_lower[i], point[i]);
            _upper[i] = std::max(_upper[i], point[i]);
        }
    }

    /// Enlarges the box by 'eps' to each side
    void extendByEps(float eps){
        for (int i=0; i<dim; i++) {
            _lower[i] -= eps;
            _upper[i] += eps;
        }
    }


private:

	C  _lower;
	C  _upper;
};

#endif // BOX_HH_
