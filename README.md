# hamplicons - estimating indel sizes in amplicons using Hamming distances

Back to length testing reads after categorizing them based on approx 30 bp match in start and end
also additional matching. very important
much faster than de2 and reporting as good or better than de3.
reason for leaving de3 with its bowtie alignment was that bowtie occasionally fails with larger indels.

Thinking about how to deal with similar PCR products, e.g. alternative primers for same product
ideal is a frameshifting alternative primer but we have historically used a few non-frameshifting ones
Currently I deal with by NOT allowing a ham score larger than or equal to the 2 products ham, which
for a non-frameshifted alternative primer usually means 1.
The issue with this is that if there naturally exists a SNP in the early part, I cant approve the product.
to get around this I currently do a ultra-relaxed matching where these will be found, and then hammed to find out which
product they are most similar to.
And alternative is to figure out the differences.
E.g. AAATGCC
VS   TAATGCC
now i know that I cant allow a SNP on pos 0, but  SNPS are ok elsewhere.
to generalize you would create a binary mask like 0111111 and then when hamming you also give the mask and no hamming is done on 0
