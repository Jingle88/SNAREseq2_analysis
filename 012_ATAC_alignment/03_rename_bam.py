import pysam
import re
from optparse import OptionParser

def add_cb_tags(input_bam, output_bam):
    """
    从QNAME中提取CB和UB信息，添加为标准标签
    """
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        with pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
            for read in infile:
                # 提取CB和UB信息
                match = re.search(r'_CB:Z:([^_]+)_UB:Z:([^_]+)$', read.query_name)
                if match:
                    cb = match.group(1)
                    ub = match.group(2)

                    # 清理QNAME
                    clean_qname = re.sub(r'_CB:[^_]+_UB:[^_]+$', '', read.query_name)
                    read.query_name = clean_qname

                    # 添加CB和UB标签
                    read.set_tag("CB", cb, "Z")
                    read.set_tag("UB", ub, "Z")

                outfile.write(read)

def main():
    parser = OptionParser(usage="usage: %prog [options] input.bam output.bam")
    parser.add_option("-i", "--index", action="store_true", dest="create_index", 
                      default=False, help="Create index for output BAM file")
    
    (options, args) = parser.parse_args()
    
    if len(args) != 2:
        parser.error("Incorrect number of arguments")
    
    input_bam, output_bam = args
    
    add_cb_tags(input_bam, output_bam)
    
    if options.create_index:
        pysam.index(output_bam)

if __name__ == "__main__":
    main()
