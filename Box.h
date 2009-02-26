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

#ifndef BOX_HH_
#define BOX_HH_


/** \brief A axis-parallel box in a Euclidean space
    \tparam C Type used for coordinates
    \tparam dim Dimension of the box
*/
template<typename C, int dim>
class Box
{
public:

	~Box()
	{}

	Box(const C& lower, const C& upper) : _lower(lower), _upper(upper)
    {}

	Box(const McVec3f& lower, const McVec3f& upper)
    {
        for (int i=0; i<3; i++) {
            _lower[i] = lower[i];
            _upper[i] = upper[i];
        }
    }

	Box(const Box& b) : _lower(b._lower), _upper(b._upper)
    {}

	bool contains(const C& c) const
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
