# Make FASTA contig IDs unique within a file.
#
# A FASTA record's ID is the first whitespace-delimited token of the header;
# everything after the first space is a free-text description that tools ignore.
# Some assemblies in circulation are headed
#
#     >SAMPLE.fasta 1
#     >SAMPLE.fasta 2
#
# which looks unique line-by-line but gives every contig the ID "SAMPLE.fasta".
# snippy rejects that reference outright:
#
#     Duplicate sequence SAMPLE.fasta in <ref>
#
# and it is silently worse elsewhere -- anything keying on contig ID collapses
# records together. 17 of the 112 B. pseudomallei assemblies in the reference
# collection are affected, including 134- and 227-contig drafts where all contigs
# share one ID.
#
# Only colliding IDs are rewritten, to "<id>_<n>". The first occurrence of an ID
# is always left alone and the description is preserved verbatim, so a file whose
# IDs are already unique is reproduced byte-for-byte. Renames are reported on
# stderr so a run's logs record exactly what was changed.
#
# Usage: awk -f dedupe_contig_ids.awk in.fasta > out.fasta

/^>/ {
    hdr = substr($0, 2)

    id = hdr
    sub(/[ \t].*$/, "", id)

    rest = hdr
    sub(/^[^ \t]*/, "", rest)

    if (id == "") { id = "contig" }

    if (id in seen) {
        seen[id]++
        newid = id "_" seen[id]
        while (newid in used) {
            seen[id]++
            newid = id "_" seen[id]
        }
        used[newid] = 1
        renamed++
        printf(">%s%s\n", newid, rest)
        printf("  renamed duplicate contig id %s -> %s\n", id, newid) > "/dev/stderr"
    } else {
        seen[id] = 0
        used[id] = 1
        printf(">%s%s\n", id, rest)
    }
    next
}

{ print }

END {
    if (renamed > 0) {
        printf("dedupe_contig_ids: rewrote %d duplicate contig id(s) in %s\n",
               renamed, FILENAME) > "/dev/stderr"
    }
}
