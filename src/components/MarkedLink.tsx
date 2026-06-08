import { useState } from 'react';
import './MarkedLink.css';

interface Props {
  href: string;
  label: string;
  note: string; // small badge revealed on hover/focus
  tone?: 'warn' | 'info'; // warn = rust (spoilers), info = teal (archive/etc.)
  italic?: boolean;
}

/* A link with a small note that reveals on hover/focus. The note is absolutely
   positioned so it never reflows the line. `label`/`note` are props (not
   children) so they survive island hydration. */
export default function MarkedLink({ href, label, note, tone = 'warn', italic }: Props) {
  const [hot, setHot] = useState(false);
  return (
    <span
      className={`ml ml--${tone}${hot ? ' is-hot' : ''}`}
      onMouseEnter={() => setHot(true)}
      onMouseLeave={() => setHot(false)}
    >
      <a
        className="ml__link"
        href={href}
        target="_blank"
        rel="noopener noreferrer"
        onFocus={() => setHot(true)}
        onBlur={() => setHot(false)}
      >
        {italic ? <em>{label}</em> : label}
      </a>
      <span className="ml__flag">{note}</span>
    </span>
  );
}
