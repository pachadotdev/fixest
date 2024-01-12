// Lag related functions

#include "05_0_misc.hpp"

[[cpp11::register]] integers cpp_lag_obs_(integers id, integers time,
                                          int nlag) {
  // in case of ties, we sum
  // must be two consecutive years
  // returns an observation nber of where to find the lagged obs
  // note that we forbid duplicate rows!!!!! => in case of duplicate we use a
  // rule of thumb

  int nobs = id.size();
  int obs;
  int id_current;
  int time_current;

  //  OLD: res(nobs, NA_INTEGER)

  // NEW: pre-allocate vector with NA integer
  writable::integers res(nobs);
  for (int i = 0; i < nobs; i++) {
    res[i] = NA_INTEGER;
  }

  int i, j, diff_time;

  if (nlag > 0) {
    i = 0;
    while (i < nobs) {
      // check_user_interrupt(); // this is (too) costly
      id_current = id[i];
      time_current = time[i];
      obs = i + 1;  // observation, R style
      j = i + 1;
      while (j < nobs) {
        diff_time = time[j] - time_current;
        if (id[j] != id_current) {
          // we start someone else
          i = j - 1;  // minus 1 because after we indent i
          break;
        } else if (diff_time > nlag) {
          // we are too far => stop
          break;
        } else if (diff_time == 0) {
          // it's the same time/id => we skip
          ++i;
        } else if (diff_time < nlag) {
          // not far enough
        } else if (diff_time == nlag) {
          // match!
          res[j] = obs;
        }

        ++j;
      }
      ++i;
    }
  } else if (nlag < 0) {
    /**************************************************************************
     * NOTA: I could have tweaked the previous if() to get rid of the condition
     *       but the code would have lost in clarity.
     *       For the lead: opposite to what is done before
     ***************************************************************************/
    int nlead = -nlag;
    i = nobs - 1;
    while (i >= 0) {
      // check_user_interrupt(); // this is (too) costly
      id_current = id[i];
      time_current = time[i];
      obs = i + 1;  // observation, R style
      j = i - 1;
      while (j >= 0) {
        diff_time = time_current - time[j];
        if (id[j] != id_current) {
          // we start someone else
          i = j + 1;  // plus 1 because after we dedent i
          break;
        } else if (diff_time > nlead) {
          // we are too far => stop
          break;
        } else if (diff_time == 0) {
          // it's the same time/id => we skip
          --i;
        } else if (diff_time < nlead) {
          // not far enough
        } else if (diff_time == nlead) {
          // match!
          res[j] = obs;
        }

        --j;
      }
      --i;
    }
  } else {
    for (int i = 0; i < nobs; ++i) {
      res[i] = i + 1;
    }
  }

  return (res);
}
