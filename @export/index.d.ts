// export R# source type define for javascript/typescript language
//
// package_source=CellRender

declare namespace CellRender {
   module _ {
      /**
      */
      function onLoad(): object;
   }
   /**
   */
   function annotation_workflow(): object;
   /**
   */
   function assemble_metabolic_graph(app: any, context: any): object;
   /**
   */
   function assemble_transcript_graph(app: any, context: any): object;
   /**
   */
   function compile_genbank(cad_registry: any, gbff: any): object;
   /**
     * @param workdir default value Is ``./``.
     * @param upstream_size default value Is ``150``.
     * @param tag_genbank_accid default value Is ``false``.
     * @param verbose default value Is ``true``.
   */
   function extract_gbff(src: any, workdir?: any, upstream_size?: any, tag_genbank_accid?: any, verbose?: any): object;
   /**
   */
   function extract_gene_table(app: any, context: any): object;
   /**
     * @param outputdir default value Is ``null``.
     * @param up_len default value Is ``150``.
     * @param biocyc default value Is ``./biocyc``.
     * @param regprecise default value Is ``./RegPrecise.Xml``.
   */
   function modelling_cellgraph(src: any, outputdir?: any, up_len?: any, biocyc?: any, regprecise?: any): object;
   /**
     * @param outputdir default value Is ``null``.
   */
   function modelling_kinetics(src: any, outputdir?: any): object;
   /**
     * @param host default value Is ``localhost``.
     * @param port default value Is ``3306``.
   */
   function open_registry(user: any, passwd: any, host?: any, port?: any): object;
   /**
   */
   function tfbs_motif_scanning(app: any, context: any): object;
}
