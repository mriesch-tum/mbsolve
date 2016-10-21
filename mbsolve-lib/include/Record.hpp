#ifndef RECORD_H
#define RECORD_H

#include <string>

namespace mbsolve {

class Record
{
public:
    /* what */
    /* string enum?? */
    /* check content where? */
    std::string name;

    /* TODO where */
    /* spatial: Quantity x sample */
    /* TODO when */
    /* temporal: Quantity dt inverval */
    /* or certain times */


    virtual bool record(unsigned int index) const = 0;


    /* virtual delgate function*/
    /* function(real *dst, real * src, unsigned int size) */
    /* save() { call function( } */

    /* then: test_and_save() { if (record) save() } */


    /* TODO: make Result 2D? */
    /* TODO: Result assignment operator???? */
    /* TODO: use library for matrices? */
};

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
