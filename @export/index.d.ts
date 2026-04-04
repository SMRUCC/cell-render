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
     * @param knn default value Is ``200``.
     * @param workdir default value Is ``./``.
   */
   function analysis(embedding_file: any, knn?: any, workdir?: any): object;
   /**
   */
   function annotation_workflow(): object;
   /**
     * @param diamond default value Is ``Call "Sys.which"("diamond")``.
     * @param n_threads default value Is ``32``.
   */
   function batch_diamond(source_dir: any, result_dir: any, diamond?: any, n_threads?: any): object;
   /**
   */
   function build_project(app: any, context: any): object;
   /**
   */
   function check_build_module(flag: any): object;
   /**
     * @param vcell_name default value Is ``null``.
   */
   function compile_model(proj_file: any, save_model: any, registry: any, vcell_name?: any): object;
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
   function extract_genomes(src: any, outputdir: any): object;
   /**
     * @param outputdir default value Is ``./``.
   */
   function genome_scatter_viz(data: any, outputdir?: any): object;
   /**
   */
   function list_batch_models(): object;
   /**
   */
   function make_blastp_term(proj_file: any, model_dir: any): object;
   /**
     * @param diamond default value Is ``Call "Sys.which"("diamond")``.
   */
   function make_diamond(local_db: any, diamond?: any): object;
   /**
   */
   function make_diamond_hits(app: any, context: any): object;
   /**
   */
   function make_genbank_proj(app: any, context: any): object;
   /**
     * @param workdir default value Is ````.
     * @param batch_process default value Is ``false``.
   */
   function make_genbank_proj_file(src: any, release_dir: any, workdir?: any, batch_process?: any): object;
   /**
   */
   function make_terms(app: any, context: any): object;
   /**
   */
   function make_TRN(app: any, context: any): object;
   /**
   */
   function model_accession_id(replicons: any): object;
   /**
     * @param outputdir default value Is ``null``.
     * @param name default value Is ``null``.
     * @param up_len default value Is ``150``.
     * @param localdb default value Is ``null``.
     * @param diamond default value Is ``Call "Sys.which"("diamond")``.
     * @param domain default value Is ``Call "c"("bacteria",
     *       "plant",
     *       "animal",
     *       "fungi")``.
     * @param builds default value Is ``Call "c"("TRN_network", "Metabolic_network")``.
     * @param n_threads default value Is ``32``.
   */
   function modelling_cellgraph(src: any, outputdir?: any, name?: any, up_len?: any, localdb?: any, diamond?: any, domain?: any, builds?: any, n_threads?: any): object;
   /**
     * @param diamond default value Is ``Call "Sys.which"("diamond")``.
     * @param n_threads default value Is ``32``.
     * @param skip_blastp default value Is ``false``.
   */
   function pangenome_analysis(src: any, result_dir: any, diamond?: any, n_threads?: any, skip_blastp?: any): object;
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
