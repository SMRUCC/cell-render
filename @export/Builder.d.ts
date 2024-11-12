﻿// export R# package module type define for javascript/typescript language
//
//    imports "Builder" from "CellRender";
//
// ref=CellRender.Builder@CellRender, Version=1.0.0.0, Culture=neutral, PublicKeyToken=null

/**
 * Helper functions for build virtualcell model
 * 
*/
declare namespace Builder {
   /**
     * @param env default value Is ``null``.
   */
   function create_modelfile(register: object, genes: any, env?: object): object;
}
