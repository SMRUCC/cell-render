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
     * @param family default value Is ``null``.
     * @param identities_cutoff default value Is ``0.8``.
     * @param minW default value Is ``0.85``.
     * @param top default value Is ``3``.
     * @param permutation default value Is ``2500``.
     * @param tqdm_bar default value Is ``true``.
     * @param env default value Is ``null``.
   */
   function motif_search(db: object, search_regions: any, family?: any, identities_cutoff?: number, minW?: number, top?: object, permutation?: object, tqdm_bar?: boolean, env?: object): object;
   /**
     * @param env default value Is ``null``.
   */
   function open_motifdb(file: any, env?: object): any;
}
