// export R# package module type define for javascript/typescript language
//
//    imports "cad_lab" from "CellRender";
//
// ref=CellRender.CadLab@CellRender, Version=1.0.0.0, Culture=neutral, PublicKeyToken=null

/**
 * 
*/
declare namespace cad_lab {
   /**
    * save the molecule expression data into the database
    *  
    *  The molecule expression data is a @``T:SMRUCC.genomics.Analysis.HTS.DataFrame.Matrix`` object, which contains the molecule expression
    *  data for each gene in the virtual cell. The @``T:CellRender.cad_lab`` object is used to access the database
    *  and save the data.
    *  
    *  The function will create a new experiment project if it does not exist, and then save the molecule
    *  expression data into the database.
    * 
    * 
     * @param cad_lab -
     * @param exp_id the experiment project id
     * @param dynaimics -
     * @param env 
     * + default value Is ``null``.
   */
   function save_expression(cad_lab: object, exp_id: string, dynaimics: object, env?: object): any;
}
