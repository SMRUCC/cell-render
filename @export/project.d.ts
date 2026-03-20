// export R# package module type define for javascript/typescript language
//
//    imports "project" from "CellRender";
//
// ref=CellRender.ProjectBuilder@CellRender, Version=1.0.0.0, Culture=neutral, PublicKeyToken=null

/**
 * 
*/
declare namespace project {
   /**
     * @param env default value Is ``null``.
   */
   function load(file: any, env?: object): any;
   /**
    * 
    * 
     * @param replicons a vector of the ncbi genbank object of the genome replicons
     * @param env -
     * 
     * + default value Is ``null``.
   */
   function new(replicons: any, env?: object): object;
   /**
     * @param env default value Is ``null``.
   */
   function save(proj: object, file: any, env?: object): any;
}
