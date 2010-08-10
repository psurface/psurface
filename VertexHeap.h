#ifndef VERTEX_HEAP
#define VERTEX_HEAP

#include <vector>

class VertexHeap{

public:

    class ErrorValue {
    public:

        float value;

        bool blocked:1;

    public:
        ErrorValue() : value(0), blocked(false){}

//      ErrorValue(float err, bool b) : featureLineStatus(0) {
//          value   = err;
//          blocked = b;
//      }



        void block() {blocked = true;}

        void unblock() {blocked = false;}

        bool isBlocked() const {
            return blocked;
        }

        friend int operator > (const ErrorValue& a, const ErrorValue& b) {
            if (a.isBlocked() == b.isBlocked())
                return a.value > b.value;
            else 
                return (a.isBlocked()) ? 1 : 0;
        }

        friend int operator < (const ErrorValue& a, const ErrorValue& b) {
            if (a.isBlocked() == b.isBlocked())
                return a.value < b.value;
            else 
                return (a.isBlocked()) ? 0 : 1;
        }
    };

    class HeapEntry {
    public:

        int vertex;
        ErrorValue    error;

        //////////////////////////////////
        
    };


public:

    void print() {
        for (int i=0; i<heapSize; i++)
            printf("Nr. %d,   value = %f,   blocked = %d\n", i, array[i].error.value, array[i].error.isBlocked());
    }

    void buildHeap(const std::vector<ErrorValue>& errors) {

        heapSize = errors.size();
        array.resize(heapSize);
        vertexNumbers.resize(heapSize);

        // enter vertices into array in any order
        int i;
        //for (cV=inList.first(); cV; cV=inList.succ(cV), i++) {
        for (i=0; i<heapSize; i++) {
            array[i].vertex = i;
            //cV->number = i;
            vertexNumbers[i] = i;
            array[i].error = errors[i];
        }

        // heapify the whole thing
        for (i= (heapSize-2)/2; i>=0; i--)
            heapify(i);
    }

    ///
    float getMinError() const {
        assert(array.size());
        return array[0].error.value;
    }

    ///
    ErrorValue getMinErrorStatus() const {
        assert(array.size());
        return array[0].error;
    }

    ///
    bool isBlockedMin() const {
        if (!heapSize)
            return true;

        return array[0].error.isBlocked();
    }
        
    ///
    int getMin() const {
        if (array.size())
            return array[0].vertex;
        else
            return -1;
    }

    ErrorValue getError(int v) const {
        return array[vertexNumbers[v]].error;
    }

    ///
    int extractMin() {

        //assert(isHeap(0));

        if (!array.size())
            return -1;

        int min = array[0].vertex;
        array[0] = array[heapSize-1];
        //array[0].vertex->number = 0;
        vertexNumbers[array[0].vertex] = 0;

        heapSize--;
        heapify(0);
        return min;
    }

    ///
    void insert(int v, const ErrorValue& err) {
        
        int i = heapSize;
        heapSize++;

        while (i>0 && array[parent(i)].error > err) {
            array[i] = array[parent(i)];
            //array[i].vertex->number = i;
            vertexNumbers[array[i].vertex] = i;

            i = parent(i);
        }

        array[i].vertex  = v;
        array[i].error   = err;
        //v->number = i;
        vertexNumbers[v] = i;
    }

    void reposition(int v, const ErrorValue& newError) {

        //assert(array[v->number].vertex==v);
        assert(array[vertexNumbers[v]].vertex == v);

        //ErrorValue oldError = array[v->number].error;
        ErrorValue oldError = array[vertexNumbers[v]].error;

        if (newError>oldError){
            //array[v->number].error   = newError;
            array[vertexNumbers[v]].error = newError;
            //heapify(v->number);
            heapify(vertexNumbers[v]);
        } else {
         
            //int i = v->number;
            int i = vertexNumbers[v];

            while (i>0 && array[parent(i)].error > newError) {
                array[i] = array[parent(i)];
                //array[i].vertex->number = i;
                vertexNumbers[array[i].vertex] = i;

                i = parent(i);
            }
        

        array[i].vertex = v;
        array[i].error  = newError;
        //v->number = i;
        vertexNumbers[v] = i;
        }
    }
        
protected:


    void heapify(int i) {

        int l = left(i);
        int r = right(i);
        int smallest;

        if (l<heapSize && array[l].error < array[i].error)
            smallest = l;
        else
            smallest = i;
        
        if (r<heapSize && array[r].error < array[smallest].error)
            smallest = r;

        if (smallest!=i){
            exchange(i, smallest);
            heapify(smallest);
        }
    }
            




    int parent(int i) const {return (i-1)/2;}

    int left(int i) const {return 2*i+1;}

    int right(int i) const {return 2*i+2;}

    void exchange(int k, int l) {
        std::swap(array[l], array[k]);

        vertexNumbers[array[k].vertex] = k;
        vertexNumbers[array[l].vertex] = l;
        
    }

    int heapSize;

    std::vector<HeapEntry> array;

    std::vector<int> vertexNumbers;

};

#endif
