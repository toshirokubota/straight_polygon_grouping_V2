#include "Array.h"

// Default constructor for class Array
template <class Item>
Array<Item> ::Array( int size )
{
#ifdef DEBUG_Array
  cout << "constructing an array of " << size << " elements" << endl;
#endif
  Size = ( size > 0 ? size : 1 ); 
  Data = new Item[ Size ]; 
  assert( Data != 0 );
}

// Copy constructor for class Array
// must receive a reference to prevent infinite recursion
template <class Item>
Array<Item>::Array( const Array<Item>& init )
{
  Size = init.GetSize();
#ifdef DEBUG_Array
  cout << "copying an array of " << Size << " elements" << endl;
#endif
  Data = new Item [ Size ]; // create space for array
  assert( Data != 0 );    // terminate if memory not allocated
  for ( int i = 0; i < Size; i++ )
    Data[ i ] = init[ i ];  // copy init into object
}

// Destructor for class Array
template <class Item>
Array<Item> ::~Array()
{
#ifdef DEBUG_Array
  cout << "destroying an array of " << Size << " elements" << endl;
#endif
  delete [] Data;            // reclaim space for array
}


// Overloaded assignment operator
// const return avoids: ( a1 = a2 ) = a3
template <class Item>
const Array<Item>& Array<Item>::operator=( const Array<Item>& right )
{
  if ( &right != this ) {  // check for self-assignment
      
    // for arrays of different sizes, deallocate original
    // left side array, then allocate new left side array.
    if ( Size != right.GetSize() ) {
      if (Size == 1)
        delete Data;
      else
        delete [] Data;         // reclaim space
      Size = right.GetSize();     // resize this object
      Data = new Item [ Size ]; // create space for array copy
      assert( Data != 0 );    // terminate if not allocated
    }
#ifdef DEBUG_Array
    cout << "assigning an array of " << Size << " elements" << endl;
#endif

    for ( int i = 0; i < Size; i++ )
      Data[ i ] = right[ i ];  // copy array into object
  }

  return *this;   // enables x = y = z;
}


// Overloaded input operator for class Array;
// inputs values for entire array.
template <class Item>
istream &operator>>( istream& input, Array<Item>& a )
{
  for ( int i = 0; i < a.GetSize(); i++ )
    input >> a[ i ];

  return input;   // enables cin >> x >> y;
}

// Overloaded output operator for class Array 
template <class Item>
ostream &operator<<( ostream& output, const Array<Item>& a )
{
  int i, size;

  size = a.GetSize();
  for ( i = 0; i < size; i++ ) {
    output << a[ i ];
    if(i < size - 1)
      output << ",";

    if ( ( i + 1 ) % Array<Item>:: ARRAY_PRINT_NUM_PER_LINE == 0 )
      output << endl;
  }

  if ( i % Array<Item>:: ARRAY_PRINT_NUM_PER_LINE != 0 )
    output << endl;

  return output;   // enables cout << x << y;
}
