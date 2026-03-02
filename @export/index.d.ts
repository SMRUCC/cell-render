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
     * @param workdir default value Is ``./``.
   */
   function analysis(embedding_file: any, workdir?: any): object;
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
     * @param union_contigs default value Is ``250``.
   */
   function diamond_embedding(diamond_result: any, workdir?: any, union_contigs?: any): object;
   /**
     * @param workdir default value Is ``./``.
     * @param node_equals default value Is ``0.999``.
   */
   function export_tree(embedding_file: any, workdir?: any, node_equals?: any): object;
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
     * @param outputdir default value Is ``./``.
   */
   function genome_scatter_viz(data: any, outputdir?: any): object;
   /**
   */
   function local_protDb(cad_registry: any, dbfile: any): object;
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
   function open_cadlab(user: any, passwd: any, host?: any, port?: any): object;
   /**
     * @param host default value Is ``localhost``.
     * @param port default value Is ``3306``.
   */
   function open_registry(user: any, passwd: any, host?: any, port?: any): object;
   /**
   */
   function reaction_pool(cad_registry: any, repo: any): object;
   /**
     * @param workdir default value Is ``./``.
   */
   function singlecells_analysis(embedding_file: any, workdir?: any): object;
   /**
     * @param outputdir default value Is ``./``.
   */
   function singlecells_viz(rawdata: any, outputdir?: any): object;
   /**
   */
   function tfbs_motif_scanning(app: any, context: any): object;
}
