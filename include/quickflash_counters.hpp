// C++ header file quickflash_counters.hpp

/*
  By Nathan C. Hearn
     October 30, 2006

  Declares classes for implementing different types of counters.

  Class VectorCounter
  ===================

  The input of the constructor specifies the number of "elements" along
  each axis.  The size of axisSize is the number of dimensions for the
  counter space.  Along each axis, the values of the counter run from zero
  to (axisSize[axis] - 1).  The counter is initially set to the zero vector.
  (In this case, the initial setting -- the zero vector -- is referred to as
  the vector startVector.)

  An additional constructor allows one to offset the starting values for
  each axis with the vector startVector.  axisSize still specifies the number
  of elements along each axis, but the values of the counter along a 
  particular axis will run from startVector[axis] to (startVector[axis] 
  + axisSize[axis] - 1).  The counter is initially set to startVector.

  Incrementing the counter proceeds by incrementing the least significant
  element of the counter (counterState[dims - 1]).  If this operation
  causes the element to reach the value of (axisSize[element] 
  + startVector[element]), where element is equal to (dims - 1), the 
  counterState[element] is set to startVector[element].  The next more-
  siginificant element is located (the value of "element" is decremented), 
  and counterState[element] is incremented.  If this element of counterState 
  reaches the value of (axisSize[element] + startVector[element]), then it
  is set to startVector[element], and then the process is repeated with the 
  more-significant elements.  If the most-significant digit goes out of 
  bounds, the counter is considered out of bounds, and can no longer be 
  incremented.

  Internally, the state of the counter is recorded in relativeCounterState.
  Each element of relativeCounterState always runs from zero to 
  (axisSize[element] - 1).  When the user asks for the state of the counter
  the object will report the quantity (relativeCounterState + startVector).
*/

#ifndef QUICKFLASH_COUNTERS_HPP
#define QUICKFLASH_COUNTERS_HPP


#include <vector>


namespace QuickFlash
{
class VectorCounter
  {
  public:
    VectorCounter();
    VectorCounter(const std::vector<unsigned int> & limits);
    VectorCounter(const VectorCounter & source);
    ~VectorCounter() { }

    VectorCounter & operator=(const VectorCounter & source)
      {
      reset(source);
      return *this;
      }

    void reset();
    void reset(const std::vector<unsigned int> & limits);
    void reset(const VectorCounter & source);

    unsigned int operator[](const int axis_index) const
      { return current_state[axis_index]; }

    void get_state(std::vector<unsigned int> & state) const
      { state = current_state; }

    const std::vector<unsigned int> & get_state() const 
      { return current_state; }

    unsigned int get_dims() const { return counter_dims.size(); }

    int increment();

    void reset_counter();

    void move_out_of_bounds()
      {
      current_state = counter_dims;
      counter_in_bounds = false;
      }

    bool in_bounds() const { return counter_in_bounds; }

    unsigned int get_num_states() const { return num_states; }

    unsigned int get_scalar_index() const;

    unsigned int get_scalar_index_reverse() const;

  private:
    std::vector<unsigned int> counter_dims;
    std::vector<unsigned int> current_state;

    unsigned int num_states;

    bool counter_in_bounds;
  };


}  // End namespace QuickFlash


#endif  // QUICKFLASH_COUNTERS_HPP
