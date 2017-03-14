#ifndef RECORD_H
#define RECORD_H

#include <string>
#include <Quantity.hpp>

namespace mbsolve {

enum RecordType { HField, EField, Density, D22, D11 };

class Record
{
private:
public:
    Record(const std::string& name,
	   const RecordType& type,
	   unsigned int i,
	   unsigned int j,
	   const Quantity& interval = 0.0,
	   const Quantity& position = -1.0) :
	Name(name), Type(type), I(i), J(j),
	Position(position), Interval(interval) { }

    std::string Name;

    unsigned int I, J;

    RecordType Type;

    Quantity Position;

    Quantity Interval;
};

/* virtual bool record(unsigned int index) const = 0;*/

/*
  RecordWrtTime


  inline record(index) return true;





  RecordWrtPosition
  size = N_x * N_t/Intervals

  bool record(index) return index % intervals == 0;


 */

/* initialize
   foreach rec in records {
   Result res(rec.x, rec.y);
   res.Record = &rec;

   if (rec.name == "")
   rec.src_addr = ...;
   else...


   results.push_back(res);

   }
 */

/* in loop:


//   foreach rec in records {
   foreach res in results {
   if res->Record.record(i) {
   determine src addr?
   copy(src, rec.result.get_row_pointer(res->Record.get_row(i)), res->Record.rowsize)

   }
   }
 */


}

#endif
