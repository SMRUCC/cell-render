// export R# package module type define for javascript/typescript language
//
//    imports "workflow" from "CellRender";
//
// ref=CellRender.workflow@CellRender, Version=1.0.0.0, Culture=neutral, PublicKeyToken=null

/**
 * annotation workflow
 * 
*/
declare namespace workflow {
   /**
   */
   function open_datapool(dir: string): object;
   /**
    * extract of the protein fasta sequence data to file
    * 
    * 
     * @param proj -
     * @param file -
     * @param env -
     * 
     * + default value Is ``null``.
   */
   function save_proteins(proj: object, file: any, env?: object): any;
   /**
     * @param env default value Is ``null``.
   */
   function set_blastp_result(blastp_hits: any, proj: object, group: string, env?: object): any;
   /**
     * @param env default value Is ``null``.
   */
   function set_tfbs(proj: object, tfbs: any, env?: object): any;
   /**
    * extract of the tss upstream location site sequence data
    * 
    * 
     * @param proj -
   */
   function tss_upstream(proj: object): object;
}
