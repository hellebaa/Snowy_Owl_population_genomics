# Putative Centromere Annotations

Add manually curated putative centromere intervals in:

```text
data/metadata/centromeres/putative_centromeres.bed
```

The analysis scripts accept either a tab-separated, headerless BED-like file or a headered centromere table.

Required columns:

```text
chromosome  start0  end
```

Optional columns:

```text
centromere_id  morphology  evidence  notes
```

Example:

```text
SUPER_3  12300000  17600000  SUPER_3_cen  metacentric  family0_plus_cytology  large internal satellite block
SUPER_20  1  2200000  SUPER_20_cen  acrocentric  family0_plus_cytology  terminal satellite block
```

Coordinates should be BED-style: 0-based start, half-open end. The matching pixy script converts pixy 1-based windows to BED-style intervals before overlap.

Alternatively, a headered table can use 1-based inclusive coordinates:

```text
Chr  Start  End  Chr_size  Morphology_new
SUPER_3  12440487  16991552  134863544  Meta
SUPER_20  13251143  17541233  17580604  Acro
```

For this headered format, `Start` is converted internally to BED-style `start0 = Start - 1`.

By default, pericentromeric regions are defined as 1 Mb upstream and downstream of each annotated centromere interval, excluding the centromere interval itself. The flank size is configured in `config/paths.yml` as `pericentromere_flank_bp`.
