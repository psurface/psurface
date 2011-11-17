/**
 * @file
 * @brief dimension independent octree functionality
 */
#ifndef MULTI_DIM_OCTREE_HH
#define MULTI_DIM_OCTREE_HH

#include <iostream>
#include <vector>
#include <deque>
#include <map>
#include "Box.h"

// can be defined if desired
#ifndef MEMINCREMENT
#define MEMINCREMENT 15
#endif

namespace psurface {

/** This class implements a dimension independent structure suitable for point
 *  location. It works like a quadtree or an octree, hence the name.
 *  Items of type T can be inserted. The octree provides a
 *  method @c lookup which returns a list of candidate cells a given point or box
 *  might be contained in.
 *  In order for items to not have to be of a particular type or to implement a
 *  particular interface the needed geometric information is provided by an
 *  associated functor object. The functor class F has to provide the operator
 *
 *  @code
 *  bool operator()(const CoordType& lower, const CoordType& upper, const T& item);
 *  @endcode
 *
 *  This method should return TRUE if the item intersects the box, FALSE otherwise.
 *
 *  Elements to be inserted into the octree are not copied. Instead, a
 *  pointer to the original element is stored. Consequently, elements must
 *  not be moved or deleted after they have been inserted into the tree.
 *
 *  While removal of items is perfectly supported, exhaustive use of this feature is
 *  not recommended. This is due to the tree refining on insert but not coarsening on
 *  removal. If many items are removed from the leaf cells branches are not collapsed
 *  even though maybe being sufficiently sparse.
 *
 *  \tparam T the data type of the stored data
 *  \tparam F the data type of the above mentioned functor; this class has to support operator() (see class description)
 *  \tparam C the type of the coordinates used; (const) access via operator[] is required
 *  \tparam dim the dimension of the tree (has to be the same as for the coordinates)
 *  \tparam uniform_access If TRUE, the user claims that access to the items geometry is uniform, and so only one
 *  functor is stored, which is, of course, much more efficient! For more information see insert
 */
template <class T, typename F, typename C, int dim, bool uniform_access = false>
class MultiDimOctree
{
private:
	typedef std::map<const T*, const F*>   FunctorMapType;


public:

	/// @brief the type of the stored data
	typedef T                 DataType;

	/// @brief the type of the functor that tells
	/// about the items' geometry
	typedef F                 CoordFunctor;

	/// @brief a type describing a box in dim dimensions
	typedef Box<C, dim>       BoxType;

	/// @brief the type of the coordinates used
	typedef C                 CoordType;

	/// @brief the type of the container that stores the
	/// results of lookup operations
	typedef std::vector<T*>   ResultContainer;

	/// @brief a constant depicting the # of subcells
	/// a cell is devided into
	static const int SUBCELLS = 1 << dim;

	/** Constructor. The argument @c bbox specifies the spatial region
	 * covered by the octree. During insertion of elements this region
	 * is successively subdivided into smaller subregions. Elements which
	 * do not intersect the original domain specified by @c bbox will not
	 * be inserted into any leaf of the octree.
   *
	 * The optional argument @c maxDepth denotes the maximum depth of
	 * the octree, while @c maxElemPerLeaf denotes the maximum number of
	 * elements stored in a leaf node. If this number is exceeded after
	 * a new element has been inserted and if the maximum depth of the
	 * octree has not yet been reached, 2^dim new leafs are created and the
	 * elements are distributed among the new leafs. A single element
	 * might be inserted into multiple leafs if the intersection test
	 * passes multiple times.
	 * @param bbox the domain of the octree
	 * @param maxDepth the maximum # of refined levels
	 * @param maxElemPerLeaf if reached in a cell the cell will be subdivided
   */
	MultiDimOctree(const BoxType &bbox, int maxDepth=6, int maxElemPerLeaf=10);

	/// @brief Default constructor.
	MultiDimOctree();

	/// @brief Destructor (frees all memory).
	virtual ~MultiDimOctree();

	/** Inserts a single element into the octree. Elements are not
   *  copied, but are referenced via a pointer. Therefore, elements must
   *  not be moved or deleted after insertion.
   *  Note that if uniform_access was specified TRUE, that each element's
   *  geometric data must be accessible in the same way, i.e. with the functor
   *  object given here. Otherwise a map is managed that holds the functor for
   *  each single item inserted. When performing lookups this obviously is a
   *  drawback since for each item touched the map has to be questioned to
   *  have the functor decide about the item being in the query box or not.
   *  @param element the to-insert element
   *  @param f_ the functor associated with the element
   */
	bool insert(T* element, F* f_);

	/**
	 * @brief removes an element from the octree
	 * @param element the element to remove
	 */
	bool remove(T* element)
	{
		if (!this->loadFunctor(element) || !remove(0, box, element))
		{
			return false;
		}
		else
		{
                    mapping.erase(element);
                    return true;
		}
	}

	/**@name lookup methods */
	//@{
	/** This method appends all elements which potentially may contain
        point @c pos to the dynamic array @c result. The array is not cleared
        in advance allowing you to collect results for multiple points. */
    int lookup(const std::tr1::array<C,dim>& pos, ResultContainer& result);

	/** Same as lookup except that indices instead of pointers are
        returned. Requires prior call to enableUniqueLookup. */
	int lookupIndex(const std::tr1::array<C,dim>& pos, std::vector<int>& result);

	/** This methods appends all elements that intersect a given box. */
	int lookup(const BoxType &queryBox, ResultContainer& result);

	/** Same as lookup except that indices instead of pointers are
        returned. Requires prior call to enableUniqueLookup. */
	int lookupIndex(const BoxType& queryBox, std::vector<int>& result);
	//@}

	/// Removes all elements and deletes all leafs of the octree.
	void clear();

	/// Calls @c clear and initializes octree from scratch.
	void init(const BoxType &bbox, int maxDepth=6, int maxElemPerLeaf=10);

	/// Print some statistics to stdout.
	void info();

	/// Returns size of complete octree in bytes.
	int memSize();

	/// Returns true if octree contains no elements.
	int isEmpty() const { return (allElements.front().n==0); }

	/// Returns maximum depth of octree.
	int getMaxDepth() const { return maxDepth; }

	/** Sets maximum depth of octree. The tree is not restructured
        in this call, so call this method before inserting elements. */
	void setMaxDepth(int val)	{ maxDepth = val; }

	/// Returns maximum number of elements per leaf.
	int getMaxElemPerLeaf() const { return maxElemPerLeaf; }

	/** Sets maximum number of elements per leaf. The tree is not restructured
        in this call, so call this method before inserting elements. */
	void setMaxElemPerLeaf(int val)	{ maxElemPerLeaf = val; }

	/** Since elements may be inserted into multiple leafs, it is possible
        that a single element is reported multiple times by lookup.
        This behavior can be suppressed if all elements inserted into
        the octree are arranged subsequently in a single array. In this
        case a bitfield of the size of the array is allocated. The bitfield
        is used to mark an element the first time it is found.

        The details: @c baseAddress denotes the address of the first element
        of the array, while @c nElements denotes the total size of the
        array. The array index of some element @c elem is computed using
        pointer arithmetic via <tt>&elem - baseAddress</tt>. */
	void enableUniqueLookup(int nElements, const T* baseAddress);

	/// This methods disables unique lookup of octree elements.
	void disableUniqueLookup();

	/// Returns current base address.
	const T* getBaseAddress() const { return baseAddress; }

	/// Returns global bounding box as defined in constructor or @c init.
	void getBoundingBox(BoxType &bb) const { bb = box; }

	/** for each cell call the virtual function workProc,
        which may be overloaded by derived classes. Returns 1 if operation was interrupted.*/
	int iterateCells(int leafsOnly=1);

	/// Called for each (leaf) cell (elem) from iterateCells. The depth
	/// and box as well as a (virtual) 3D cell index on this depth is
	/// provided. returns 1 to break iteration
	virtual int workProc(int elem, int depth, const BoxType &elemBox, int indices[dim])	{ return 0; }

protected:

	int iterateCells(int elem, int depth, const BoxType &elemBox,
			int leafsOnly, int indices[dim]);

	/*
	 * No description necessary since protected anyway.
	 * Short: Type of a cell that - if a leaf - stores an
	 * array of inserted items und - if not a leaf - only
	 * stores an index that helps identifying the subcells
	 * in the global array of all cells.
	 */
	struct Element
	{
		unsigned int isLeaf:1;
		// for leaves: #data items in array "indices"
		// for branches: index of first subcell in global "allElement" array
		unsigned int n:31;
		T** indices;

		Element() : isLeaf(1), n(0), indices(NULL)
		{}

		~Element()
		{
			if (indices)
				free(indices);
		}

		void remove(int remN, std::vector<bool> &remElems)
		{
			int oldN = n;
			n -= remN;
			if ((n % MEMINCREMENT) == 0)
			{ // decrease size of memory
				T** oldIndices = indices;
				indices = (T**) malloc(n*sizeof(T*));
				for (int oldI=0, i=0; oldI<oldN; oldI++)
				{
					if (!remElems[oldI])
					{
						indices[i] = oldIndices[oldI]; i++;
					}
				}
				free(oldIndices);
			}
			else { // just move items (same size of memory)
				for (int oldI=0, i=0; oldI<oldN; oldI++)
				{
					if (!remElems[oldI])
					{
						indices[i] = indices[oldI]; i++;
					}
				}
			}
		}
	};

	/// @brief random access container storing all cells in the octree
	std::deque<Element> allElements;

	bool insert(int elem, int depth, const BoxType &elemBox, T* idx);

	bool remove(int elem, const BoxType &elemBox, const T* toBeDeleted);

	void lookup(int elem, BoxType &elemBox, const std::tr1::array<C,dim>& pos, ResultContainer& result);

	void lookup(int elem, const BoxType &elemBox, const BoxType& queryBox, ResultContainer& result);

	void subdivide(int elem, const BoxType &elemBox);

	/**
	 * fills this->f with the functor for the given item
	 * @param item the item
	 * @return TRUE if uniform_access is TRUE or if the map contains
	 * a functor associated with the item
	 */
	bool loadFunctor(const T* item)
	{
		if (uniform_access)
			return this->f != NULL;
		else
		{
			typename FunctorMapType::iterator it = this->mapping.find(item);
			if (it != this->mapping.end())
			{
				this->f = (*it).second;
				return true;
			}
			else
				return false;
		}
	}


	BoxType                      box;
	int                          maxDepth;
	unsigned int                 maxElemPerLeaf;
	const T*                     baseAddress;
	std::vector<bool>            lookupFlags;
	// only needed if !uniform_access:
	// this map stores the associated functors for the items
	FunctorMapType               mapping;
	// currently "loaded" functor
	const F*                     f;
};

/// @if EXCLUDETHIS

template <class T, typename F, typename C, int dim, bool uniform_access>
MultiDimOctree<T, F, C, dim, uniform_access>::MultiDimOctree(const BoxType &box, int maxDepth, int maxElemPerLeaf) : box(box)
{
	init(box, maxDepth, maxElemPerLeaf);
}


template <class T, typename F, typename C, int dim, bool uniform_access>
MultiDimOctree<T, F, C, dim, uniform_access>::MultiDimOctree()
{
	baseAddress = 0;
	maxDepth = 0;
	maxElemPerLeaf = 0;
	allElements.clear();
	allElements.push_back(Element());
}


template <class T, typename F, typename C, int dim, bool uniform_access>
MultiDimOctree<T, F, C, dim, uniform_access>::~MultiDimOctree()
{
	// Empty, new std::vector frees all memory automatically
}


template <class T, typename F, typename C, int dim, bool uniform_access>
void MultiDimOctree<T, F, C, dim, uniform_access>::enableUniqueLookup(int n, const T* addr)
{
	baseAddress = addr;
	lookupFlags.resize(n);
	for (size_t i=0; i<lookupFlags.size(); i++)
		lookupFlags[i] = false;
}


template <class T, typename F, typename C, int dim, bool uniform_access>
void MultiDimOctree<T, F, C, dim, uniform_access>::disableUniqueLookup()
{
	baseAddress = 0;
	lookupFlags.resize(0);
}


template <class T, typename F, typename C, int dim, bool uniform_access>
void MultiDimOctree<T, F, C, dim, uniform_access>::clear()
{
	baseAddress = 0;
	lookupFlags.resize(0);

	// Just keep the root node. remax(1,1) is wrong since the root node
	// then would not be marked as a leaf.
	allElements.clear();
	allElements.push_back(Element());
}


template <class T, typename F, typename C, int dim, bool uniform_access>
void MultiDimOctree<T, F, C, dim, uniform_access>::init(const BoxType& b, int depth, int elemPerLeaf)
{
	box = b;

	baseAddress = 0;
	lookupFlags.resize(0);
	maxDepth = depth;
	maxElemPerLeaf = elemPerLeaf;

	// Create root node.
	allElements.clear();
	allElements.push_back(Element());
}


template <class T, typename F, typename C, int dim, bool uniform_access>
bool MultiDimOctree<T, F, C, dim, uniform_access>::insert(T* element, F* f_)
{
	if (uniform_access)
	{
		// set the functor
		this->f = f_;
		if (this->f != NULL && (*this->f)(this->box.lower(), this->box.upper(), *element))
		{
			// return the success of the insert method
			return insert(0, 0 , box, element);
		}
	}
	else
	{
		// check if the element is already inserted and whether it lies
		// inside the octree's bounding box
		if (f_ != NULL && this->mapping.find(element) == this->mapping.end() &&
				(*f_)(this->box.lower(), this->box.upper(), *element))
		{
			this->mapping[element] = f_;
			if (!insert(0, 0 , box, element))
			{
				this->mapping.erase(element);
				return false;
			}
			else
				return true;
		}
		else
			return false;
	}
	return false;
}





template <class T, typename F, typename C, int dim, bool uniform_access>
bool MultiDimOctree<T, F, C, dim, uniform_access>::insert(int elem, int depth, const BoxType &elemBox, T* idx)
{
	// if element is leaf -> simply insert and be done!
	Element& element = allElements[elem];
	if (element.isLeaf)
	{

		// but if element already contains max. number of items
		// subdivide and put items in subcells...
		// ...and then insert the new item (goto DESCEND)
		if (depth<maxDepth && element.n>=maxElemPerLeaf)
		{
			subdivide(elem, elemBox);
			goto DESCEND;
		}

		if (element.n % MEMINCREMENT == 0)
		{
			int newSize = element.n + MEMINCREMENT;
			if (element.indices)
			{
				element.indices =	(T**) realloc(element.indices, newSize*sizeof(T*));
			}
			else
			{
				element.indices = (T**) malloc(newSize*sizeof(T*));
			}
		}

		element.indices[element.n++] = idx;
		return true;
	}

	// this element is a parent, so go visit its children
	// and insert there
	DESCEND:

	// get the functor of the element
	if (!this->loadFunctor(idx))
		return false; // if not found in map => abort!

	int firstChild = element.n;

	depth++;

	// the return value
	bool inserted = false;

	// helpful points that describe box corners
        std::tr1::array<C,dim> upper, lower;
	// iterate over all subcells and check for intersections between
	// item and subcells
	for (int j = 0; j < SUBCELLS; ++j)
	{
		// compute the next child cell
		for(int i = 0; i < dim; ++i)
		{
			if (j & (1 << i))
			{
				lower[i] = elemBox.center()[i];
				upper[i] = elemBox.upper()[i];
			}
			else
			{
				lower[i] = elemBox.lower()[i];
				upper[i] = elemBox.center()[i];
			}
		}
		// if the item intersects the child element's box
		// insert it into this box
		BoxType childElemBox(lower, upper);
		if ((*f)(lower, upper, *idx))
			inserted = inserted || insert(firstChild+j, depth, childElemBox, idx);
	}
	return inserted;
}


template <class T, typename F, typename C, int dim, bool uniform_access>
int MultiDimOctree<T, F, C, dim, uniform_access>::iterateCells(int leafsOnly)
{
	int indices[dim];
	for (int i = 0; i < dim; ++i)
		indices[i] = 0;
	return iterateCells(0, 0, box, leafsOnly, indices);
}


template <class T, typename F, typename C, int dim, bool uniform_access>
int MultiDimOctree<T, F, C, dim, uniform_access>::iterateCells(int elem, int depth, const BoxType &elemBox,
		int leafsOnly, int indices[dim])
{
	Element& element = allElements[elem];
	if (element.isLeaf || !leafsOnly)
	{
		if (workProc(elem, depth, elemBox, indices))
			return 1;
	}
	if (element.isLeaf)
		return 0;

	int firstChild = element.n;

	std::tr1::array<C,dim> lower, upper, center = elemBox.center();

	// prepare all indices for the next stage
	depth++;
	int temp_indices[dim];
	for (int i = 0; i < dim; ++i)
		indices[i] *= 2;

	for (int j = 0; j < SUBCELLS; ++j)
	{
		// reset the temporary indices
		for (int i = 0; i < dim; ++i)
			temp_indices[i] = indices[i];

		// compute the next child cell
		for(int i = 0; i < dim; ++i)
		{
			if (j & (1 << i))
			{
				lower[i] = elemBox.center()[i];
				upper[i] = elemBox.upper()[i];
				temp_indices[i] += 1;
			}
			else
			{
				lower[i] = elemBox.lower()[i];
				upper[i] = elemBox.center()[i];
			}
		}

		BoxType childElemBox(lower, upper);
		if (iterateCells(firstChild+j, depth, childElemBox, leafsOnly, temp_indices))
			return 1;
	}
	return 0;
}


template <class T, typename F, typename C, int dim, bool uniform_access>
void MultiDimOctree<T, F, C, dim, uniform_access>::subdivide(int elem, const BoxType &elemBox)
{
	int childIndex = allElements.size();

	Element& element = allElements[elem];

	int nIndices = element.n;

	element.isLeaf = 0;
	element.n = childIndex;

	for (int i = 0; i < SUBCELLS; i++)
	{
		allElements.push_back(Element());
		// should not be necessary because these are default
		// vaules set by the constructor anyway
		//allElements.back().isLeaf = 1;
		//allElements.back().n = 0;
		//allElements.back().indices = 0;
	}

	// insert again into current element but since this is
	// not a leaf anymore the elements are put into the
	// newly created subcells properly
	for (int i = 0; i < nIndices; i++)
		insert(elem, 999, elemBox, element.indices[i]);

	// not a leaf anymore, so no items stored here either!
	if (element.indices)
	{
//		delete[] element.indices;
		free(element.indices);
		element.indices = 0;
	}
}


template <class T, typename F, typename C, int dim, bool uniform_access>
int MultiDimOctree<T, F, C, dim, uniform_access>::memSize()
{
	int bytes = allElements.size() * sizeof(Element);
	for (typename std::deque<Element>::reverse_iterator rit = allElements.rbegin(); rit != allElements.rend(); ++rit)
	{
		Element& element = *rit;
		if (element.isLeaf)
		{
			int n = element.n / MEMINCREMENT;
			bytes += (n+1)*MEMINCREMENT*sizeof(int);
		}
	}
	return bytes;
}


template <class T, typename F, typename C, int dim, bool uniform_access>
void MultiDimOctree<T, F, C, dim, uniform_access>::info()
{
	int nNodes = allElements.size();
	int nLeafs = 0;
	int nElements = 0;
	int minNumElements = 999999;
	int maxNumElements = 0;

	for (typename std::deque<Element>::iterator it = allElements.begin(); it != allElements.end(); ++it)
	{
		Element& node = *it;
		if (node.isLeaf)
		{
			nLeafs++;
			int n = node.n;
			if (n < minNumElements)
				minNumElements = n;
			if (n > maxNumElements)
				maxNumElements = n;
			nElements += n;
		}
	}

	std::cout << "MultiDimOctree: " << nNodes-nLeafs << " nodes,"
	          << nLeafs << " leafs (" << (float) memSize()/(1024*1024) << " MB)" << std::endl;
	std::cout << "MultiDimOctree: " << nElements << " elements,"
	          << " (" << (float)nElements/nLeafs << " per leaf, " << minNumElements << "..." << maxNumElements << ")" << std::endl;
}


template <class T, typename F, typename C, int dim, bool uniform_access>
int MultiDimOctree<T, F, C, dim, uniform_access>::lookup(const std::tr1::array<C,dim>& pos, ResultContainer& result)
{
	BoxType b(box);

	if (b.contains(pos))
		lookup(0, b, pos, result);

	// Ensure that lookupFlags is clean again...
	if (baseAddress)
	{
		for (int i=result.size()-1; i>=0; i--)
		{
			int n = result[i] - baseAddress;
			lookupFlags[n] = false;
		}
	}
	return result.size();
}


template <class T, typename F, typename C, int dim, bool uniform_access>
int MultiDimOctree<T, F, C, dim, uniform_access>::lookup(const BoxType& queryBox, ResultContainer& result)
{
	BoxType b(box);

	if (b.intersects(queryBox))
		lookup(0, b, queryBox, result);

	// Ensure that lookupFlags is clean again...
	if (baseAddress)
	{
		for (int i=result.size()-1; i>=0; i--)
		{
			int n = result[i] - baseAddress;
			lookupFlags[n] = false;
		}
	}
	return result.size();
}


template <class T, typename F, typename C, int dim, bool uniform_access>
int MultiDimOctree<T, F, C, dim, uniform_access>::lookupIndex(const BoxType& queryBox, std::vector<int>& result)
{
	ResultContainer tmpResult;
	lookup(queryBox, tmpResult);

	int n = tmpResult.size();
	for (int i=0; i<n; i++)
		result.push_back(tmpResult[i] - baseAddress);
	return result.size();
}


template <class T, typename F, typename C, int dim, bool uniform_access>
void MultiDimOctree<T, F, C, dim, uniform_access>::lookup(int elem, const BoxType &elemBox, const BoxType& queryBox, ResultContainer& result)
{
	Element& element = allElements[elem];

	if (element.isLeaf)
	{
		for (unsigned int i=0; i<element.n; i++)
		{
			T* t = element.indices[i];

			// get the functor of the element and check for intersection
			if (this->loadFunctor(t) && (*this->f)(queryBox.lower(), queryBox.upper(), *t))
			{
				if (baseAddress)
				{ // this indicates unique lookup strategy
					int k = t - baseAddress;
					if (lookupFlags[k] == 0)
					{
						result.push_back(t);
						lookupFlags[k] = true;
					}
				} else
					result.push_back(t); // t may be appended multiple times
			}
		}
	}
	else
	{
		int firstChild = element.n;

		// These intersection tests are faster than the generic BoxType member function,
		// since we do not have to check lower and upper coordinates of both boxes

		// the boundary of the next subcell is stored in here
                std::tr1::array<C,dim> lower, upper;
		for (int j = 0; j < SUBCELLS; ++j)
		{
			bool intersects_subcell = true;
			// compute the next child cell
			for(int i = 0; i < dim; ++i)
			{
				if (j & (1 << i))
				{
					lower[i] = elemBox.center()[i];
					upper[i] = elemBox.upper()[i];
					intersects_subcell = intersects_subcell && (queryBox.upper()[i] >= lower[i]);
				}
				else
				{
					lower[i] = elemBox.lower()[i];
					upper[i] = elemBox.center()[i];
					intersects_subcell = intersects_subcell && (queryBox.lower()[i] < upper[i]);
				}
			}
			// if intersecting then descend recursively to subcell
			if (intersects_subcell)
			{
				BoxType childElemBox(lower, upper);
				lookup(firstChild+j, childElemBox, queryBox, result);
			}
		}
	}
}


template <class T, typename F, typename C, int dim, bool uniform_access>
int MultiDimOctree<T, F, C, dim, uniform_access>::lookupIndex(const std::tr1::array<C,dim>& pos, std::vector<int>& result)
{
	ResultContainer tmpResult;
	lookup(pos, tmpResult);

	int n = tmpResult.size();
	for (int i=0; i<n; i++)
		result.push_back(tmpResult[i] - baseAddress);
	return result.size();
}


template <class T, typename F, typename C, int dim, bool uniform_access>
void MultiDimOctree<T, F, C, dim, uniform_access>::lookup(int elem, BoxType &elemBox, const std::tr1::array<C,dim>& pos, ResultContainer& result)
{
	Element& element = allElements[elem];

	if (element.isLeaf)
	{
		for (unsigned int i=0; i<element.n; i++)
		{
			T* t = element.indices[i];
			if (baseAddress)
			{ // this indicates unique lookup strategy
				unsigned int k = t - baseAddress;
				if (lookupFlags[k] == 0)
				{
					result.push_back(t);
					lookupFlags[k] = true;
				}
			}
			else
				result.push_back(t); // t may be appended multiple times
		}
	}
	else
	{
		int firstChild = element.n;

		std::tr1::array<C,dim> lower, upper;

		int config = 0;
		// compute the subcell in which the point is located
		for (int i = 0; i < dim; ++i)
		{
			if (pos[i] >= elemBox.center()[i])
			{
				lower[i] = elemBox.center()[i];
				upper[i] = elemBox.upper()[i];
				config |= (1 << i);
			}
			else
			{
				lower[i] = elemBox.lower()[i];
				upper[i] = elemBox.center()[i];
			}
		}
		BoxType childElemBox(lower, upper);
		lookup(firstChild+config, childElemBox, pos, result);
	}
}


template <class T, typename F, typename C, int dim, bool uniform_access>
bool MultiDimOctree<T, F, C, dim, uniform_access>::remove(int elem, const BoxType &elemBox, const T* toBeDeleted)
{
	Element& element = allElements[elem];

	if (element.isLeaf)
	{
		std::vector<bool> remElems(element.n, false);
		const T* t;
		int remN = 0;
		// always remove all appearances of elemPtr
		for (unsigned int i=0; i<element.n; i++)
		{
			t = element.indices[i];
			if (t == toBeDeleted)
			{
				remElems[i] = true; remN++;
			}
		}
		if (remN)
		{
			element.remove(remN, remElems);
			return true;
		}
		return false;
	}
	else
	{
		int firstChild = element.n;

		// the result value
		bool removed = false;

		// helpful points that describe box corners
		std::tr1::array<C,dim> upper, lower;
		// iterate over all subcells and check for intersections between
		// item and subcells
		for (int j = 0; j < SUBCELLS; ++j)
		{
			// compute the next child cell
			for(int i = 0; i < dim; ++i)
			{
				if (j & (1 << i))
				{
					lower[i] = elemBox.center()[i];
					upper[i] = elemBox.upper()[i];
				}
				else
				{
					lower[i] = elemBox.lower()[i];
					upper[i] = elemBox.center()[i];
				}
			}
			// if the item intersects the child element's box
			// insert it into this box
			BoxType childElemBox(lower, upper);
			if ((*this->f)(lower, upper, *toBeDeleted))
				removed = removed || remove(firstChild+j, childElemBox, toBeDeleted);
		}

		return removed;
	}
}

/// @endif

} // namespace psurface

#endif // MULTI_DIM_OCTREE_HH

/// @}
