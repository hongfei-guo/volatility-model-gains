# Source data

The third-party index data are not included in this archive. They were
obtained from Yahoo Finance and remain subject to the provider's terms of use
and any restrictions imposed by its data suppliers. Replicators must obtain
the source observations directly from Yahoo Finance or another licensed
provider.

Yahoo's current terms and historical-data documentation are available at:

- <https://legal.yahoo.com/us/en/yahoo/terms/otos/index.html>
- <https://in.help.yahoo.com/kb/SLN2311.html>

## Required files

Create `data/source/` and place one daily CSV file for each market in it:

| File | Index | Yahoo ticker | First required output date | Last date | Output observations |
|---|---|---|---:|---:|---:|
| `DAX.csv` | DAX Performance Index | `^GDAXI` | 2000-01-04 | 2025-12-30 | 6,600 |
| `FTSE.csv` | FTSE 100 | `^FTSE` | 2000-01-05 | 2025-12-31 | 6,566 |
| `NIKKEI.csv` | Nikkei 225 | `^N225` | 2000-01-05 | 2025-12-30 | 6,368 |
| `SHCOMP.csv` | SSE Composite | `000001.SS` | 2000-01-05 | 2025-12-31 | 6,294 |
| `SP500.csv` | S&P 500 | `^GSPC` | 2000-01-04 | 2025-12-31 | 6,538 |

Each file must contain `Date` and `Close` columns; capitalization does not
matter. Include at least one trading day before the first required output date
so that its return can be calculated. A common request window beginning in
December 1999 is sufficient.

## Return construction

Run the following command from the archive root:

```text
Rscript code/build_returns.R
```

The script sorts each market's observations by date and calculates daily
percentage log returns as

`100 × [log(close_t) − log(close_{t−1})]`.

It checks the required market calendars and writes the analysis series to
`reproduced/returns/`. Missing observations are not interpolated or
forward-filled. Provider revisions can cause newly obtained historical data to
differ from the observations used in the paper.
