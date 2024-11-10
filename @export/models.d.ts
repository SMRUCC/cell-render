// export R# package module type define for javascript/typescript language
//
//    imports "models" from "biocad_registry";
//
// ref=biocad_registry.models@biocad_registry, Version=1.0.0.0, Culture=neutral, PublicKeyToken=null

/**
*/
declare namespace models {
   /**
   */
   function enzyme(biocad_registry: object, ec_number: string): any;
   /**
     * @param env default value Is ``null``.
   */
   function enzyme_function(biocad_registry: object, enzyme_id: string, ec_number: string, metabolites: any, env?: object): any;
   /**
   */
   function subcellular_location(biocad_registry: object, name: string, topology: string): object;
   /**
     * @param desc default value Is ``''``.
   */
   function vocabulary_id(biocad_registry: object, term: string, category: string, desc?: string): object;
}
