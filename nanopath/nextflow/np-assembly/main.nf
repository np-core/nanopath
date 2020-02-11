#!/usr/bin/env nextflow

/*

Pipeline            np-assembly
Version             v0.1

Description         Nanopore workflow for long-read genome assembly and polishing

Authors             Eike Steinig

*/

log.info """
--------------------------------------

outdir          : $params.outdir

ont             : $params.ont
quality         : $params.quality
length          : $params.length
subsample       : $params.subsample
assembler       : $params.assembler
genome_size     : $params.genome_size
opts            : $params.opts
medaka_model    : $params.medaka_model
illumina        : $params.illumina
assembler       : $params.illumina_assembler
--------------------------------------
"""


fastq_nanopore = Channel
    .fromPath(params.ont)
    .map { file -> tuple(file.baseName, file) }

process NanoQC {
    
    tag { id }
    label "ont"

    publishDir "$params.outdir/fastq", mode: "copy", pattern: "*.txt"

    input:
    set id, file(fq) from fastq_nanopore

    output:
    set id, file("${id}.filtered.fq") into rasusa_fastq
    set file("${id}.filtered.stats.txt"), file("${id}.prefiltered.stats.txt")

    """
    NanoStat --threads $task.cpus --fastq $fq > ${id}.prefiltered.stats.txt
    NanoFilt --quality $params.quality --length $params.length $fq > ${id}.filtered.fq
    NanoStat --threads $task.cpus --fastq ${id}.filtered.fq > ${id}.filtered.stats.txt
    """

}

process Rasusa {
    
    tag { id }
    label "ont"

    publishDir "$params.outdir/fastq", mode: "copy"

    input:
    set id, file(fq) from rasusa_fastq

    output:
    set id, file("${id}.sub.fq") into assembly_fastq
    
    script:

    if ( params.subsample > 0 )
        
        """
        rasusa -c $params.subsample -g $params.genome_size -i $fq > ${id}.sub.fq
        """
    else 
        // sym link here?

        """
        cp $fq ${id}.sub.fq
        """

}

process Assembly {
    
  tag { "$id:$params.assembler" }
  label "assembly"

  publishDir "$params.outdir/assembly", mode: "copy"

  input:
  set id, file(fq) from assembly_fastq

  output:
  file("${id}.assembly_info.txt")
  set id, file("${id}.assembly.fasta"), file(fq) into (racon_assembly, pilon_assembly)
  set id, file("${id}.assembly_graph.gfa") into bandage_assembly
  
  script:

  if ( params.assembler == 'flye' )

    """
    flye --nano-raw $fq --genome-size $params.genome_size $params.opts -t $task.cpus -o assembly
    mv assembly/assembly_info.txt ${id}.assembly_info.txt
    mv assembly/assembly.fasta ${id}.assembly.fasta
    mv assembly/assembly_graph.gfa ${id}.assembly_graph.gfa
    """

}

process Racon {

    tag { id }
    label "racon"

    input:
    set id, file(assembly), file(fastq) from racon_assembly

    output:
    set id, file("${id}.racon.fasta"), file(fastq) into medaka_racon

    script:
    """
    minimap2 -x map-ont -t $task.cpus $assembly $fastq > assembly_1.paf
    racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_1.paf $assembly > assembly_consensus_1.fasta
    minimap2 -x map-ont -t $task.cpus assembly_consensus_1.fasta $fastq > assembly_2.paf
    racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_2.paf assembly_consensus_1.fasta > assembly_consensus_2.fasta
    minimap2 -x map-ont -t $task.cpus assembly_consensus_2.fasta $fastq > assembly_3.paf
    racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_3.paf assembly_consensus_2.fasta > assembly_consensus_3.fasta
    minimap2 -x map-ont -t $task.cpus assembly_consensus_3.fasta $fastq > assembly_4.paf
    racon -m 8 -x -6 -g -8 -w 500 -t $task.cpus $fastq assembly_4.paf assembly_consensus_3.fasta > ${id}.racon.fasta
    """
}

process Medaka {

    tag { id }
    label "medaka"

    publishDir "$params.outdir/medaka", mode: "copy"

    input:
    set id, file(racon_assembly), file(fastq) from medaka_racon

    output:
    set id, file("${id}.consensus.fasta") into medaka_consensus

    """
    medaka_consensus -i $fastq -d $racon_assembly -o racon_medaka -t $task.cpus -m $params.medaka_model
    mv racon_medaka/consensus.fasta ./${id}.consensus.fasta
    """

}


// Short read section

if (params.illumina) {

    // Pair assemblies with parsed Illumina files by ID:


    illumina_reads = Channel.fromFilePairs(params.illumina, flat: true)

    process Trimmomatic {

        tag { rid }
        label "trimmomatic"

        input:
        set rid, file(forward), file(reverse) from illumina_reads

        output:
        set rid, file("${rid}_1P.fq.gz"), file("${rid}_2P.fq.gz") into (pilon_illumina, shovill_illumina, pilon_assembly_illumina)

        """
        trimmomatic PE $forward $reverse \
        -threads $task.cpus -phred33 -baseout ${rid}.fq.gz \
        ILLUMINACLIP:$baseDir/resources/trimmomatic/all_adapters.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

    }
    
    medaka_pilon = medaka_consensus.cross(pilon_illumina)
                    .map { crossed ->
                        if (crossed[0][0] == crossed[1][0]){
                            tuple( crossed[0][0], crossed[0][1],crossed[1][1], crossed[1][2] )
                        } else {
                            null
                        }
                    }
                    .filter { it != null }
    
    process MedakaPilon {

        tag { aid }
        label "pilon"

        publishDir "$params.outdir/pilon_medaka", mode: "copy"

        input:
        set aid, file(consensus_assembly), file(forward), file(reverse) from medaka_pilon

        output:
        file("${aid}.medaka.pilon.fasta")


        """
        minimap2 -ax sr $consensus_assembly $forward $reverse > aln.sam && \
            samtools view -S -b aln.sam > aln.bam && samtools sort aln.bam -o alignment1.bam && \
            samtools index alignment1.bam

        pilon --genome $consensus_assembly --frags alignment1.bam --outdir correction1 --changes
        
        minimap2 -ax sr correction1/pilon.fasta $forward $reverse > aln.sam && \
            samtools view -S -b aln.sam > aln.bam && samtools sort aln.bam -o alignment2.bam && \
            samtools index alignment2.bam

        pilon --genome correction1/pilon.fasta  --frags alignment2.bam --outdir correction2 --changes && \
        mv correction2/pilon.fasta ${aid}.medaka.pilon.fasta
        """

    }

    process Shovill {

        tag { rid }
        label "shovill"

        publishDir "$params.outdir/$params.illumina_assembler", mode: "copy"

        input:
        set rid, file(forward), file(reverse) from shovill_illumina

        output:
        file("${rid}.fasta")

        """
        shovill --R1 ${forward} --R2 ${reverse} --cpus $task.cpus --ram $task.memory \
        --depth 100 --assembler $params.illumina_assembler --outdir $rid --force
        mv ${rid}/contigs.fa ${rid}.fasta
        """

    }   
    
    assembly_pilon = pilon_assembly.cross(pilon_assembly_illumina)
                    .map { crossed ->
                        if (crossed[0][0] == crossed[1][0]){
                            tuple( crossed[0][0], crossed[0][1], crossed[1][1], crossed[1][2] )
                        } else {
                            null
                        }
                    }
                    .filter { it != null }
    
    process AssemblyPilon {

        tag { aid }
        label "pilon"

        publishDir "$params.outdir/pilon_assembly", mode: "copy"

        input:
        set aid, file(assembly), file(forward), file(reverse) from assembly_pilon

        output:
        file("${aid}.assembly.pilon.fasta")


        """
        minimap2 -ax sr $assembly $forward $reverse > aln.sam && \
            samtools view -S -b aln.sam > aln.bam && samtools sort aln.bam -o alignment1.bam && \
            samtools index alignment1.bam

        pilon --genome $assembly --frags alignment1.bam --outdir correction1 --changes
        
        minimap2 -ax sr correction1/pilon.fasta $forward $reverse > aln.sam && \
            samtools view -S -b aln.sam > aln.bam && samtools sort aln.bam -o alignment2.bam && \
            samtools index alignment2.bam

        pilon --genome correction1/pilon.fasta  --frags alignment2.bam --outdir correction2 --changes && \
        mv correction2/pilon.fasta ${aid}.assembly.pilon.fasta
        """

    }

}