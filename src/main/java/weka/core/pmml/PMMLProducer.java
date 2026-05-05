/*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */

/*
 *    PMMLProducer.java
 *    Copyright (C) 2014 University of Waikato, Hamilton, New Zealand
 *
 */

package weka.core.pmml;

import weka.core.Instances;

/**
 * Interface to something that can produce a PMML representation of itself.
 * 
 * @author David Persons
 * @author Mark Hall (mhall{[at]}pentaho{[dot]}com)
 * @version $Revision: $
 */
public interface PMMLProducer {

  /**
   * Produce a PMML representation
   * 
   * @param train the training data that might have been used by the
   *          implementer. If it is not needed by the implementer then clients
   *          can safely pass in null
   * 
   * @return a string containing the PMML representation
   */
  String toPMML(Instances train);
}
